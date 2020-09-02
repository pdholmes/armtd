classdef test_preslice_rotatotope_RTD_planner < robot_rotatotope_RTD_planner
    %ROBOT_ROTATOTOPE_RTD_PLANNER general class for using rotatotope
    % representations of the forward reachable set to plan safe robot
    % trajectories
    
    properties
    end
    
    methods
        function P = test_preslice_rotatotope_RTD_planner(varargin)
            P@robot_rotatotope_RTD_planner('fetch_rotatotope_FRS', varargin{:});
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
            
            % generate FRS
            %             P.R = eval(sprintf('%s(q_0, q_dot_0, P.FRS_options)', P.rotatotope_FRS_class));
            
            %%% TESTING: try fixing joints 3 - 5 at constant velocity, only
            %%% optimize over joints 1, 2, and 6.
            pre_slice_dim = {[4]; [4]; [4;3]; [4;3]; [4;3]; [4]};  % for each JRS, dimensions to pre-slice
            pre_slice_values = {[q_dot_0(1)]; [q_dot_0(2)]; [q_dot_0(3); 0]; [q_dot_0(4); 0]; [q_dot_0(5); 0]; [q_dot_0(6)]}; % for each JRS, values to pre-slice
            k_dim = {[3]; [3]; []; []; []; [3]}; % for each JRS, dimensions of trajectory parameters
            k_names = {{'ka1'}; {'ka2'}; {}; {}; {}; {'ka6'}}; % store list of all parameters
            P.R = fetch_rotatotope_FRS(q_0, q_dot_0, P.FRS_options, ...
                'pre_slice_dim', pre_slice_dim, ...
                'pre_slice_values', pre_slice_values, ...
                'k_dim', k_dim, ...
                'k_names', k_names);
            
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
                %%% MALICIOUS HACK!!
                grad_h_joint = [-grad_q_min([1, 2, 6], :), grad_q_max([1, 2, 6], :), -grad_q_dot_min([1, 2, 6], :), grad_q_dot_max([1, 2, 6], :)];
                gradh = [gradh, grad_h_joint];
            end
                        
            %%% Obstacle constraint generation:
            for i = 1:length(P.R.A_con) % for each obstacle
                for j = 1:length(P.R.A_con{i}) % for each link
                    idx = find(~cellfun('isempty', P.R.A_con{i}{j}));
                    for k = 1:length(idx) % for each time step
%                         n_k = size(P.R.k_con{i}{j}{idx(k)}, 1);
                        k_idx = [];
                        for l = 1:length(P.R.k_con_list{i}{j}{idx(k)})
                            k_idx(l, 1) = find(strcmp(P.R.k_con_list{i}{j}{idx(k)}{l}, P.R.k_list));
                        end
                        [h_obs, grad_h_obs_tmp] = P.R.evaluate_constraint(P.R.A_con{i}{j}{idx(k)}, P.R.b_con{i}{j}{idx(k)}, P.R.k_con{i}{j}{idx(k)}, P.R.k_con_list{i}{j}{idx(k)}, k_opt(k_idx, 1));
                        h = [h; h_obs + epsilon]; % add epsilon to h_obs to turn <= into strictly <
                        grad_h_obs = zeros(length(k_opt), size(grad_h_obs_tmp, 2));
                        grad_h_obs(k_idx, :) = grad_h_obs_tmp;
                        if nargout > 2
                            gradh = [gradh, grad_h_obs];
                        end
                    end
                end
            end
            
            %%% Self-intersection constraint generation:
            for i = 1:length(P.R.A_con_self) % for each pair of joints that can intersect
                idx = find(~cellfun('isempty', P.R.A_con_self{i}));
                for j = 1:length(idx) % for each (nonempty) time step
%                     n_k = size(P.R.k_con_self{i}{idx(j)}, 1);
                    k_idx = [];
                    for l = 1:length(P.R.k_con_self_list{i}{idx(j)})
                        k_idx(l, 1) = find(strcmp(P.R.k_con_self_list{i}{idx(j)}{l}, P.R.k_list));
                    end
                    [h_self, grad_h_self_tmp] = P.R.evaluate_constraint(P.R.A_con_self{i}{idx(j)}, P.R.b_con_self{i}{idx(j)}, P.R.k_con_self{i}{idx(j)}, P.R.k_con_self_list{i}{idx(j)}, k_opt(k_idx, 1));
                    h = [h; h_self + epsilon]; % add epsilon to h_obs to turn <= into strictly <
                    grad_h_self = zeros(length(k_opt), size(grad_h_self_tmp, 2));
                    grad_h_self(k_idx, :) = grad_h_self_tmp;
                    if nargout > 2
                        gradh = [gradh, grad_h_self];
                    end
                end
            end
        end
        
        function [T, U, Z] = generate_trajectory(P, q_0, q_dot_0, k_opt)
            
            %%% MALICIOUS HACK:
            k_opt = [k_opt(1); k_opt(2); 0; 0; 0; k_opt(3)];
            
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
            
            %%% MALICIOUS HACK:
            k = [k(1); k(2); 0; 0; 0; k(3)];
            q_plan = q_0 + q_dot_0*P.t_move + (1/2)*k*P.t_move^2;
        end
        
        function [q_min, q_max, q_dot_min, q_dot_max, grad_q_min, grad_q_max, grad_q_dot_min, grad_q_dot_max] = compute_max_min_states(P, q_0, q_dot_0, k)
            % this function is specific to the particular trajectory
            % paramterization used in the ARMTD paper.
            
            %%% MALICIOUS HACK:
            k = [k(1); k(2); 0; 0; 0; k(3)];
            
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

