classdef robot_rotatotope_RTD_planner < robot_arm_generic_planner
    %ROBOT_ROTATOTOPE_RTD_PLANNER general class for using rotatotope
    % representations of the forward reachable set to plan safe robot
    % trajectories
    
    properties
        t_total = 1;
        time_discretization = 0.01;
        
        trajopt_start_tic;
        
        FRS_options;
        
        % for current FRS
        rotatotope_FRS_class;
        R ; % reachable set contained as robot_rotatotope_FRS class
        O ; % hold on to obstacles
        
        % from previous time step
        Z_prev = []; % previous state trajectory
        R_prev ;
        q_0_prev = [];
        q_dot_0_prev = [];
        k_opt_prev = []; 
        
        iter = 0;
        first_iter_pause_flag = true ;
        
        % for cost function
        use_end_effector_for_cost_flag = false ;
        forward_kinematics_end_effector = [] ;
        
        % for cuda
        use_cuda_flag = false;
    end
    
    methods
        function P = robot_rotatotope_RTD_planner(varargin)
            t_move = 0.5;
            lookahead_distance = 0.4;
            HLP = robot_arm_straight_line_HLP( );
            P@robot_arm_generic_planner('lookahead_distance', lookahead_distance, 't_move', t_move, 'HLP', HLP, ...
                varargin{3:end}) ;
            P.rotatotope_FRS_class = varargin{1};
            P.FRS_options = varargin{2};

            % init info object
            P.init_info()
        end
        
        function init_info(P)
            P.info = struct('T',[],'U',[],'Z',[],'waypoint',[],...
                'obstacles',[],'q_0',[],'q_dot_0',[],'k_opt',[]) ;
        end
        
        function [T,U,Z] = replan(P,agent_info,world_info)
            P.vdisp('Replanning!',5)
            
            P.vdisp('Seting joint state and speed limits',7)
            % set any joint limits that are +Inf to 200pi and -Inf to -200pi
            joint_limit_infs = isinf(P.arm_joint_state_limits) ;
            P.arm_joint_state_limits(1,joint_limit_infs(1,:)) = -200*pi ;
            P.arm_joint_state_limits(2,joint_limit_infs(2,:)) = +200*pi ;
            
            speed_limit_infs = isinf(P.arm_joint_speed_limits) ;
            P.arm_joint_speed_limits(1,speed_limit_infs(1,:)) = -200*pi ;
            P.arm_joint_speed_limits(2,speed_limit_infs(2,:)) = +200*pi ;
            
            P.vdisp('Generating cost function',6)
            
            % get forward kinematics function
            if P.use_end_effector_for_cost_flag
                P.forward_kinematics_end_effector = agent_info.get_end_effector_location ;
            end
            
            % generate a waypoint in configuration space
            if P.first_iter_pause_flag && P.iter == 0
               pause; 
            end
            P.iter = P.iter + 1;
            planning_time = tic;
            q_des = P.HLP.get_waypoint(agent_info,world_info,P.lookahead_distance) ;
            if isempty(q_des)
                P.vdisp('Waypoint creation failed! Using global goal instead.', 3)
                q_des = P.HLP.goal ;
            end
                        
            P.vdisp('Generating constraints',6)
            % get current obstacles
            P.O = world_info.obstacles;
            
            % get current state of robot
            q_0 = agent_info.state(P.arm_joint_state_indices, end) ;
            q_dot_0 = agent_info.state(P.arm_joint_speed_indices, end) ;
            
            if ~P.use_cuda_flag
                % generate FRS
                P.R = eval(sprintf('%s(q_0, q_dot_0, P.FRS_options)', P.rotatotope_FRS_class));
                
                % map obstacles to trajectory parameter space
                P.R = P.R.generate_constraints_from_obstacles(P.O);
                
                % protect against self-intersections:
                P.R = P.R.generate_self_intersection_constraints;
                
                P.vdisp('Replan is calling trajopt!',8)
                [k_opt, trajopt_failed] = P.trajopt(q_0, q_dot_0, q_des);
                P.vdisp('Processing trajopt result',8)
                
                % plotting
%                 P.R.plot();
%                 P.R.plot_slice(k_opt);
%                 pause;
            else
                P.vdisp('Replan is calling trajopt!',8)
                R_cuda = robot_arm_FRS_rotatotope_fetch_cuda(q_0, q_dot_0, q_des, P.O, zeros(6,1), P.FRS_options);
                disp(R_cuda.mex_res);
                k_opt = R_cuda.mex_res;
                if(length(R_cuda.mex_res) == 6)
                    trajopt_failed = false;
                else
                    trajopt_failed = true;
                end
                P.vdisp('Processing trajopt result',8)
            end
            
            if ~trajopt_failed
                P.vdisp('New trajectory found!',3);
                
                [T, U, Z] = P.generate_trajectory(q_0, q_dot_0, k_opt);
                
                toc(planning_time);
             
                P.R_prev = P.R;
                P.Z_prev = Z;
                
                P.q_0_prev = q_0;
                P.q_dot_0_prev = q_dot_0;
                P.k_opt_prev = k_opt;
            else % if no safe trajectory parameter found:
                % generate a braking trajectory
                P.vdisp('Unable to find new trajectory!',3)
                T = 0:P.time_discretization:P.t_total ;
                T_plan = T(1:floor(length(T)/2));
                N_T = length(T) ;
                U = zeros(P.arm_n_inputs,N_T) ;
                if ~isnan(P.q_0_prev) % if prev trajectory was not a braking traj
                    Z_brake = P.Z_prev(:, floor(length(T)/2)+1:end);
                    z = [Z_brake(1:2:end, end)'; zeros(length(q_0), 1)'];
                    Z = [Z_brake, ones(size(Z_brake, 1), length(T_plan)).*z(:)];
                    P.q_0_prev = nan;
                    P.q_dot_0_prev = nan;
                    P.k_opt_prev = nan;
                else
                    % command to stay in place
                    z = [q_0'; zeros(length(q_0), 1)'];
                    Z = ones(length(q_0) + length(q_dot_0), length(T)).*z(:);
                    P.q_0_prev = nan;
                    P.q_dot_0_prev = nan;
                    P.k_opt_prev = nan;
                end
                toc(planning_time);
            end
            
            % save info
            P.info.T = [P.info.T, {T}] ;
            P.info.U = [P.info.U, {U}] ;
            P.info.Z = [P.info.Z, {Z}] ;
            P.info.waypoint = [P.info.waypoint, {q_des}] ;
            P.info.obstacles = [P.info.obstacles, {world_info.obstacles}] ;
            P.info.q_0 = [P.info.q_0, {q_0}] ;
            P.info.q_dot_0 = [P.info.q_dot_0, {q_dot_0}] ;
            P.info.k_opt = [P.info.k_opt, {k_opt}] ;
        end
        
        function [k_opt, trajopt_failed] = trajopt(P, q_0, q_dot_0, q_des)
            % use fmincon to optimize the cost subject to constraints
            P.vdisp('Running trajopt', 3)
            
            P.trajopt_start_tic = tic ;
            
            if P.use_end_effector_for_cost_flag
                ee_des = P.forward_kinematics_end_effector(q_des) ;
                cost_func = @(k) P.eval_cost_end_effector(k, q_0, q_dot_0, ee_des);
            else
                cost_func = @(k) P.eval_cost(k, q_0, q_dot_0, q_des);
            end
            constraint_func = @(k) P.eval_constraint(k, q_0, q_dot_0);
            
            % generate upper and lower bounds
            lb = P.R.c_k - P.R.delta_k;
            ub = P.R.c_k + P.R.delta_k;
            
%             initial_guess = (lb + ub)/2;
            initial_guess = rand_range(lb, ub);
           
%             options = optimoptions('fmincon','SpecifyConstraintGradient',true);
            options = optimoptions('fmincon','SpecifyConstraintGradient',true, 'CheckGradients', true);
            [k_opt, ~, exitflag, ~] = fmincon(cost_func, initial_guess, [], [], [], [], lb, ub, constraint_func, options) ;
            
            trajopt_failed = exitflag <= 0 ;
        end
        
        function [cost] = eval_cost(P, k, q_0, q_dot_0, q_des)
            % generate a simple cost function
           q_plan = compute_q_plan(P, q_0, q_dot_0, k);
           cost = sum((q_plan - q_des).^2);
        end
        
        function c = eval_cost_end_effector(P, k, q_0, q_dot_0, ee_des)
            q_plan = compute_q_plan(P, q_0, q_dot_0, k);
            ee_plan = P.forward_kinematics_end_effector(q_plan(:,end)) ;
            c = sum((ee_plan - ee_des).^2) ;    
        end
        
        function [h, heq, gradh, gradheq] = eval_constraint(P, k_opt, q_0, q_dot_0)
            epsilon = 1e-3;
            heq = [];            
            h = [];
            if nargout > 2
                gradheq = [];
                gradh = [];
            end
            %%% Joint limit constraints:
            [q_min, q_max, q_dot_min, q_dot_max, grad_q_min, grad_q_max, grad_q_dot_min, grad_q_dot_max] = P.compute_max_min_states(q_0, q_dot_0, k_opt);

            h_joint = [];
            h_joint = [h_joint; P.arm_joint_state_limits(1, :)' - q_min];
            h_joint = [h_joint; -P.arm_joint_state_limits(2, :)' + q_max];
            h_joint = [h_joint; P.arm_joint_speed_limits(1, :)' - q_dot_min];
            h_joint = [h_joint; -P.arm_joint_speed_limits(2, :)' + q_dot_max];
                        
            h = [h; h_joint];
            if nargout > 2
                grad_h_joint = [-grad_q_min, grad_q_max, -grad_q_dot_min, grad_q_dot_max];
                gradh = [gradh, grad_h_joint];
            end
                        
            %%% Obstacle constraint generation:
            for i = 1:length(P.R.A_con) % for each obstacle
                for j = 1:length(P.R.A_con{i}) % for each link
                    idx = find(~cellfun('isempty', P.R.A_con{i}{j}));
                    for k = 1:length(idx) % for each time step
                        n_k = size(P.R.k_con{i}{j}{idx(k)}, 1);
                        [h_obs, grad_h_obs] = P.R.evaluate_constraint(P.R.A_con{i}{j}{idx(k)}, P.R.b_con{i}{j}{idx(k)}, P.R.k_con{i}{j}{idx(k)}, k_opt(1:n_k));
                        h = [h; h_obs + epsilon]; % add epsilon to h_obs to turn <= into strictly <
                        if nargout > 2
                            gradh = [gradh, [grad_h_obs; zeros(length(k_opt) - n_k, 1)]];
                        end
                    end
                end
            end
            
            %%% Self-intersection constraint generation:
            for i = 1:length(P.R.A_con_self) % for each pair of joints that can intersect
                idx = find(~cellfun('isempty', P.R.A_con_self{i}));
                for j = 1:length(idx) % for each (nonempty) time step
                    n_k = size(P.R.k_con_self{i}{idx(j)}, 1);
                    [h_self, grad_h_self] = P.R.evaluate_constraint(P.R.A_con_self{i}{idx(j)}, P.R.b_con_self{i}{idx(j)}, P.R.k_con_self{i}{idx(j)}, k_opt(1:n_k));
                    h = [h; h_self + epsilon]; % add epsilon to h_obs to turn <= into strictly <
                    if nargout > 2
                        gradh = [gradh, [grad_h_self; zeros(length(k_opt) - n_k, 1)]];
                    end
                end
            end
        end
        
        function [T, U, Z] = generate_trajectory(P, q_0, q_dot_0, k_opt)
            T = 0:P.time_discretization:P.t_total ;
            T_plan = T(1:floor(length(T)/2)+1) ;
            T_brake = T(floor(length(T)/2)+1:end) ;
            
            t_to_stop = T_brake(end) - T_brake(1);
            U = zeros(P.arm_n_inputs,length(T));
            
            q_to_peak = q_0 + q_dot_0.*T_plan + (1/2)*k_opt.*T_plan.^2;
            q_dot_to_peak = q_dot_0 + k_opt.*T_plan;
            
            q_peak = q_to_peak(:, end);
            q_dot_peak = q_dot_to_peak(:, end);
            
            T_brake = T_brake - T_brake(1);
            q_to_stop = q_peak + q_dot_peak.*T_brake + (1/2)*((0 - q_dot_peak)./t_to_stop).*T_brake.^2;
            q_dot_to_stop = q_dot_peak + ((0 - q_dot_peak)./t_to_stop).*T_brake;
            
            % remove overlap
            q_to_stop(:, 1) = [];
            q_dot_to_stop(:, 1) = [];
            
            Z = zeros(length(q_0)*2, length(T));
            for i = 1:size(q_to_peak, 1)
               Z(2*i-1:2*i, :) =  [q_to_peak(i, :), q_to_stop(i, :); q_dot_to_peak(i, :), q_dot_to_stop(i, :)];
            end
        end
        
        function [q_plan] = compute_q_plan(P, q_0, q_dot_0, k)
            % returns the configuration at t_move given initial state and
            % chosen acceleration k; note if P.t_plan = P.t_move then
            % real-time planning is enforced (more or less)
            
            q_plan = q_0 + q_dot_0*P.t_move + (1/2)*k*P.t_move^2;
        end
        
        function [q_min, q_max, q_dot_min, q_dot_max, grad_q_min, grad_q_max, grad_q_dot_min, grad_q_dot_max] = compute_max_min_states(P, q_0, q_dot_0, k)
            % this function is specific to the particular trajectory
            % paramterization used in the ARMTD paper.
            
            % compute the max and min joint positions and velocities over
            % the time horizon, given k.
            n_q = length(q_0);
            n_k = length(k);
            t_to_stop = P.t_total - P.t_move;
            
            q_peak = q_0 + q_dot_0*P.t_move + (1/2)*k*P.t_move^2;
            q_dot_peak = q_dot_0 + k*P.t_move;
            q_ddot_to_stop = ((0 - q_dot_peak)./t_to_stop);
            q_stop = q_peak + q_dot_peak.*t_to_stop + (1/2)*q_ddot_to_stop.*t_to_stop.^2;
            
            % note that position trajectories are piecewise quadratic in
            % time, and velocity trajectoreis are piecewise linear.
            % also, note that when braking, we have 0 velocity at P.t_total.
            % therefore, max and min q_dot occur at the endpoints of
            % the braking segments, or the endpoints of the "to peak"
            % segments.
            % max and min q can occur at these endpoints, or also at some
            % local max/min within the "to peak" segments.
            
            t_max_min_to_peak = -q_dot_0./k; % time of max or min for to peak dynamics
            for i = 1:n_q
                if t_max_min_to_peak(i) > 0 && t_max_min_to_peak(i) < P.t_move % this implies local max/min of q
                    q_max_min_to_peak(i, 1) = q_0(i) + q_dot_0(i)*t_max_min_to_peak(i) + (1/2)*k(i)*t_max_min_to_peak(i)^2;
                else
                    q_max_min_to_peak(i, 1) = nan;
                end
            end
            
            possible_max_min_q = [q_0, q_max_min_to_peak, q_peak, q_stop]';
            possible_max_min_q_dot = [q_dot_0, q_dot_peak, zeros(size(q_dot_0))]';
            
            [q_min, q_min_idx] = min(possible_max_min_q);
            [q_max, q_max_idx] = max(possible_max_min_q);
            [q_dot_min, q_dot_min_idx] = min(possible_max_min_q_dot);
            [q_dot_max, q_dot_max_idx] = max(possible_max_min_q_dot);
            
            % the gradients depend on where the maxes and mins occur
            possible_grad_q = [zeros(n_q, 1), ((1/2)*q_dot_0.^2)./(k.^2), ones(n_q, 1)*(1/2)*P.t_move^2, ones(n_q, 1)*((1/2)*P.t_move^2 + (1/2)*P.t_move*t_to_stop)]';
            possible_grad_q_dot = [zeros(n_q, 1), ones(n_q, 1)*P.t_move, zeros(n_q, 1)]';
            
            for i = 1:n_q
                grad_q_min(i, 1) = possible_grad_q(q_min_idx(i), i);
                grad_q_max(i, 1) = possible_grad_q(q_max_idx(i), i);
                grad_q_dot_min(i, 1) = possible_grad_q_dot(q_dot_min_idx(i), i);
                grad_q_dot_max(i, 1) = possible_grad_q_dot(q_dot_max_idx(i), i);
            end
            
            % reshape
            q_min = q_min';
            q_max = q_max';
            q_dot_min = q_dot_min';
            q_dot_max = q_dot_max';
            
            % finally, make diagonal matrices out of gradients:
            grad_q_min = diag(grad_q_min);
            grad_q_max = diag(grad_q_max);
            grad_q_dot_min = diag(grad_q_dot_min);
            grad_q_dot_max = diag(grad_q_dot_max);
        end
        
    end
end

