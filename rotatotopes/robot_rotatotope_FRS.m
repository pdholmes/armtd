classdef robot_rotatotope_FRS
    %ROBOT_ROTATOTOPE_FRS Represents the forward reachable set (FRS) of a
    %robot's links and joints as a set of rotatotopes.
    %   As described in Lemma 7, we overapproximate the FRS of a robot in
    %   workspace using sets of rotatotopes. We use these FRS's to
    %   construct collision-avoidance and self-intersection constraints,
    %   described in Sec. 4.B
    %
    %   This class should serve as a superclass, where subclasses contain
    %   robot-specific parameters.
    
    properties
        % still to do: infer these properties from a particular agent
        % class, so that robot-specific subclasses of this class are 
        % unnecessary:
        rotation_axes = []; % array of column vectors representing the body-fixed 3D rotation axis of each joint
        n_links = []; % number of physical robot links
        n_base = []; % number of "base links", not used for collision constraints but affect FRS construction
        link_predecessor_joints = {}; % cell array of vectors, each specifying all predecessor joints of each link
        link_zonotopes = {}; % zonotopes representing the volume of each link
        joint_zonotopes = {}; % zonotopes representing volume of joints at end of each link
        base_predecessor_joints = {}; % joints that precede base "link"... base link not used for collision detection, but affects FRS construction
        base_joint_zonotopes = {}; % zonotope representing volume of base joint at end of base link, used for FRS construction
        link_self_intersection = {}; % cell array of links that could intersect for self-collision constraints
        n_time_steps; % number of time steps defining the JRSs
        
        JRS_path = 'rotatotopes/joint_reachable_sets/'; % where the JRSs are stored
        JRS_key = []; % used to identify which JRS to load
        
        q_0; % initial configuration at start of planning iteration
        q_dot_0; % initial velocity at start of planning iteration
        
        % store information about trajectory parameter space:
        k_dim = []; % 
        c_k = [];
        delta_k = [];
        
        link_rotatotopes = {}; % store rotatotopes representing FRS of each link
        
        FRS_options = {}; % additional options for FRS creation
        
        % cells for storing obstacle collision-avoidance constraints
        A_con = {};
        b_con = {};
        k_con = {};
        
        % cells for storing self-intersection collision-avoidance
        % constraints
        A_con_self = {};
        b_con_self = {};
        k_con_self = {};
        
        % patch data for plotting FRSs
        plot_patch = {};
        plot_patch_exists = false;
    end
    
    methods
        function obj = robot_rotatotope_FRS(q_0, q_dot_0, FRS_options)
            %ROBOT_ROTATOTOPE_FRS Construct robot FRS represented as
            %rotatotopes.
            %   The constructor requires the initial configuration and
            %   velocity at the start of the planning iteration (q_0 and
            %   q_dot_0). Additional options can be specified in
            %   FRS_options. The constructor creates the rotatotope
            %   reachable sets of the robot links.
            obj.q_0 = q_0;
            obj.q_dot_0 = q_dot_0;
            
            JRS_key_tmp = load([obj.JRS_path, 'c_kvi.mat']);
            obj.JRS_key = JRS_key_tmp.c_kvi;
            
            if ~exist('FRS_options', 'var')
                obj.FRS_options.combs = generate_combinations_upto(200);
                obj.FRS_options.maxcombs = 200;
                obj.FRS_options.buffer_dist = 0;
                obj.FRS_options.origin_shift = zeros(3, 1);
            else
                obj.FRS_options = FRS_options;
            end
                       
        end
        
        function obj = create_FRS(obj)
            % this function compose the rotatotope forward reachable set of
            % the robot, as in Algorithm 2.
            JRS = cell(length(obj.q_dot_0), 1);
            for i = 1:length(obj.q_dot_0)
                [~, closest_idx] = min(abs(obj.q_dot_0(i) - obj.JRS_key));
                JRS_filename = sprintf('%sJRS_%0.3f.mat', obj.JRS_path, obj.JRS_key(closest_idx));
                JRS_load_tmp = load(JRS_filename);
                JRS_load{i} = JRS_load_tmp.JRS;
                
                % implementation of lines 4 and 5 of Alg. 2:
                % use A to rotate the cos and sin dimensions of JRS
                % then, slice at correct k^v_i value
                A = [cos(obj.q_0(i)), -sin(obj.q_0(i)), 0, 0, 0; sin(obj.q_0(i)), cos(obj.q_0(i)), 0, 0, 0;...
                    0, 0, 1, 0, 0; 0, 0, 0, 1, 0; 0, 0, 0, 0, 1];
                for j = 1:length(JRS_load{i})
                    JRS{j}{i} = A*zonotope_slice(JRS_load{i}{j}{1}, 4, obj.q_dot_0(i));
                end
                
                % extract K information
                obj.c_k(i, 1) = JRS{1}{i}.Z(obj.k_dim, 1);
                G = JRS{1}{i}.Z(:, 2:end);
                k_col = find(G(obj.k_dim, :) ~= 0);
                if length(k_col) == 1
                    obj.delta_k(i, 1) = G(obj.k_dim, k_col);
                else
                    error('More than one generator for k-dimension');
                end
            end
            
            obj.n_time_steps = length(JRS);
            
            % base joint rotatotope
            for i = 1:obj.n_base
                for j = 1:obj.n_time_steps
                    base_joint_rotatotopes{i}{j} = rotatotope_v2(obj.rotation_axes(:, obj.base_predecessor_joints{i}), JRS{j}(obj.base_predecessor_joints{i}), obj.base_joint_zonotopes{i});
                end
            end
            
            % link volume rotatotopes:
            for i = 1:obj.n_links
                for j = 1:obj.n_time_steps
                    link_volume_rotatotopes{i}{j} = rotatotope_v2(obj.rotation_axes(:, obj.link_predecessor_joints{i}), JRS{j}(obj.link_predecessor_joints{i}), obj.link_zonotopes{i});
                end
            end
            
            % link predecessor joint rotatotopes
            for i = 1:obj.n_links - 1
                for j = 1:obj.n_time_steps
                    link_joint_rotatotopes{i}{j} = rotatotope_v2(obj.rotation_axes(:, obj.link_predecessor_joints{i}), JRS{j}(obj.link_predecessor_joints{i}), obj.joint_zonotopes{i});
                end
            end
            
            % stack:
            % (lines 12 - 14 of Alg. 2)
            obj.link_rotatotopes = link_volume_rotatotopes;
            for i = 1:obj.n_links
                for j = 1:obj.n_time_steps
                    % stack on base
                    for k = 1:obj.n_base
                       obj.link_rotatotopes{i}{j} = obj.link_rotatotopes{i}{j}.stack(base_joint_rotatotopes{k}{j}); 
                    end
                    % stack on prev. links
                    for k = 1:i-1
                       obj.link_rotatotopes{i}{j} = obj.link_rotatotopes{i}{j}.stack(link_joint_rotatotopes{k}{j});
                    end
                end
            end
            
            % shift origins:
            % for when there is a mismatch between world frame origin and
            % base of kinematic chain for rotatotopes
            for i = 1:obj.n_links
                for j = 1:obj.n_time_steps
                    obj.link_rotatotopes{i}{j}.Vit(:, 1) = obj.link_rotatotopes{i}{j}.Vit(:, 1) + obj.FRS_options.origin_shift;
                end
            end
            
            % the following is useful for discarding unnecessary
            % constraints online. basically, we store all the vertices
            % of the trajectory parameter space, and later use them to test
            % if any choice of trajectory parameter could cause a
            % collision. if not, we will discard the constraint.
            for i = 1:obj.n_links
                cubeV = cubelet_ND(length(obj.link_predecessor_joints{i}))';
                obj.FRS_options.k_vertices{i} = (cubeV.*obj.delta_k(obj.link_predecessor_joints{i})) + obj.c_k(obj.link_predecessor_joints{i});
                obj.FRS_options.k_vertices_kappa{i} = cubeV;
            end
        end
        
        function [patch_data] = plot(obj, plot_times, colors, patch_data)
            % plots the (buffered) link_rotatotopes:
            if ~exist('colors', 'var')
                colors = repmat(1/256*[158,202,225], obj.n_links, 1);
            end
            if ~exist('plot_times', 'var')
                plot_times = 1:10:obj.n_time_steps;
            end
            if ~exist('patch_data', 'var')
                patch_data = {};
            end
            if isempty(patch_data)
                for i = 1:obj.n_links
                    for j = plot_times
                        patch_data{i}{j} = plot(obj.link_rotatotopes{i}{j}, colors(i, :), obj.FRS_options.buffer_dist);
                    end
                end
            else
                for i = 1:obj.n_links
                    for j = plot_times
                        patch_data{i}{j} = plot(obj.link_rotatotopes{i}{j}, colors(i, :), obj.FRS_options.buffer_dist, patch_data{i}{j});
                    end
                end
            end
        end
        
        function [patch_data] = plot_slice(obj, k, plot_times, colors, patch_data)
            % plots the (buffered) link_rotatotopes sliced at k:
            if ~exist('colors', 'var')
                colors = repmat([0, 1, 0], obj.n_links, 1);
            end
            if ~exist('plot_times', 'var')
                plot_times = 1:10:obj.n_time_steps;
            end
            if ~exist('patch_data', 'var')
                patch_data = {};
            end
            if isempty(patch_data)
                for i = 1:obj.n_links
                    for j = plot_times
                        patch_data{i}{j} = plot_slice(obj.link_rotatotopes{i}{j}, k(obj.link_predecessor_joints{i}), colors(i, :), obj.FRS_options.buffer_dist);
                    end
                end
            else
                for i = 1:obj.n_links
                    for j = plot_times
                        patch_data{i}{j} = plot_slice(obj.link_rotatotopes{i}{j}, k(obj.link_predecessor_joints{i}), colors(i, :), obj.FRS_options.buffer_dist, patch_data{i}{j});
                    end
                end
            end
        end
        
        function obj = generate_constraints_from_obstacles(obj, obstacles)
           for i = 1:length(obstacles)
               % get zonotope representation of obstacle:
               obstacle = obstacles{i}.zono.Z;
               obstacle_center = obstacle(:, 1);
               obstacle_generators = obstacle(:, 2:end);
               obstacle_generators(:, ~any(obstacle_generators)) = []; % delete zero columns
               for j = 1:length(obj.link_rotatotopes)
                   for k = 1:length(obj.link_rotatotopes{j})
                       % get rotatotope center and generators:
                       Vit = obj.link_rotatotopes{j}{k}.Vit;
                       Vit_center = Vit(:, 1);
                       Vit_generators = Vit(:, 2:end);
                       
                       % get k_slc and fully_slc from rotatotope
                       k_slc = obj.link_rotatotopes{j}{k}.k_slc;
                       fully_slc = obj.link_rotatotopes{j}{k}.fully_slc;
                       fully_slc_bool = (fully_slc == 1);
                       
                       % separate rotatotope generators into
                       % fully-k-sliceable vs not-fully-k-sliceable 
                       % see Lemma 9 and eqn. (19) of paper
                       Vit_slc_generators = Vit_generators(:, fully_slc_bool);
                       Vit_buf_generators = Vit_generators(:, ~fully_slc_bool);
                       
                       % create buffered obstacle (see eqn. (21))
                       O_buf = [obstacle_center - Vit_center, obstacle_generators, Vit_buf_generators, obj.FRS_options.buffer_dist/2*eye(size(Vit_center, 1))];
                       % and turn into halfpsace representation
                       [A_obs, b_obs] = polytope_PH(O_buf, obj.FRS_options);
                       
                       % used to implement eqn. (22) of paper:
                       % note that we multiply A_obs by Vit_slc... then,
                       % multiplying A_con by the corresponding generator
                       % coefficients (kappas) will yield the expression
                       % A_obs * eval(Vit_slc, k^a) found in eqn. (22)
                       A_con = A_obs*Vit_slc_generators; 
                       b_con = b_obs;
                       k_con = k_slc(:, fully_slc_bool);
                       
                       % add a test here that throws out unnecessary constraints.
                       % (there's certainly a more elegant way to do this, but we're
                       % testing all possible vertices of the parameter set K, which
                       % is a hypercube)
                       intersection_possible = 0;
                       for l = 1:size(obj.FRS_options.k_vertices{j}, 2)
                           h_obs = obj.evaluate_constraint(A_con, b_con, k_con, obj.FRS_options.k_vertices{j}(:, l));
                           if h_obs >= 0
                               intersection_possible = 1;
                               break;
                           end
                       end
                       if ~intersection_possible
                           A_con = []; b_con = []; k_con = [];
                       end
                       
                       obj.A_con{i}{j}{k} = A_con;
                       obj.b_con{i}{j}{k} = b_con;
                       obj.k_con{i}{j}{k} = k_con;
                   end
               end
           end
        end
        
        function obj = generate_self_intersection_constraints(obj)
           for i = 1:length(obj.link_self_intersection) % for each pair of links that could intersect
               % get zonotope representation of obstacle:
               for j = 1:length(obj.link_rotatotopes{1})
                   % get rotatotope center and generators:
                   Vit_1 = obj.link_rotatotopes{obj.link_self_intersection{i}(1)}{j}.Vit;
                   Vit_1_center = Vit_1(:, 1);
                   Vit_1_generators = Vit_1(:, 2:end);
                   
                   Vit_2 = obj.link_rotatotopes{obj.link_self_intersection{i}(2)}{j}.Vit;
                   Vit_2_center = Vit_2(:, 1);
                   Vit_2_generators = Vit_2(:, 2:end);

                   % get k_slc and fully_slc from rotatotope
                   k_slc_1 = obj.link_rotatotopes{obj.link_self_intersection{i}(1)}{j}.k_slc;
                   fully_slc_1 = obj.link_rotatotopes{obj.link_self_intersection{i}(1)}{j}.fully_slc;
                   fully_slc_bool_1 = (fully_slc_1 == 1);
                   
                   k_slc_2 = obj.link_rotatotopes{obj.link_self_intersection{i}(2)}{j}.k_slc;
                   fully_slc_2 = obj.link_rotatotopes{obj.link_self_intersection{i}(2)}{j}.fully_slc;
                   fully_slc_bool_2 = (fully_slc_2 == 1);

                   % separate rotatotope generators into
                   % fully-k-sliceable vs not-fully-k-sliceable 
                   % see Lemma 9 and eqn. (19) of paper
                   Vit_1_slc_generators = Vit_1_generators(:, fully_slc_bool_1);
                   Vit_1_buf_generators = Vit_1_generators(:, ~fully_slc_bool_1);
                   
                   Vit_2_slc_generators = Vit_2_generators(:, fully_slc_bool_2);
                   Vit_2_buf_generators = Vit_2_generators(:, ~fully_slc_bool_2);

                   % create buffered "link volume" (see eqns. (25) and (26))
                   % in the appendix
                   V_self = [-Vit_1_slc_generators, Vit_2_slc_generators];
                   V_buf = [Vit_1_center - Vit_2_center, Vit_1_buf_generators, Vit_2_buf_generators, 2*obj.FRS_options.buffer_dist/2*eye(size(Vit_1_center, 1))];
                   % and turn into halfpsace representation
                   [A_self, b_self] = polytope_PH(V_buf, obj.FRS_options);

                   % used to implement eqn. (27) of paper:
                   % note that we multiply A_self by Vit_self... then,
                   % multiplying A_con by the corresponding generator
                   % coefficients (kappas) will yield the expression
                   % A_self * eval(Vit_self, k^a) found in eqn. (27)
                   A_con = A_self*V_self; 
                   b_con = b_self;
                   k_con = [[k_slc_1(:, fully_slc_bool_1); nan(size(k_slc_2, 1) - size(k_slc_1, 1), size(k_slc_1(:, fully_slc_bool_1), 2))], k_slc_2(:, fully_slc_bool_2)];

                   % add a test here that throws out unnecessary constraints.
                   % (there's certainly a more elegant way to do this, but we're
                   % testing all possible vertices of the parameter set K, which
                   % is a hypercube)
                   intersection_possible = 0;
                   for l = 1:size(obj.FRS_options.k_vertices{obj.link_self_intersection{i}(2)}, 2)
                       h_obs = obj.evaluate_constraint(A_con, b_con, k_con, obj.FRS_options.k_vertices{obj.link_self_intersection{i}(2)}(:, l));
                       if h_obs >= 0
                           intersection_possible = 1;
                           break;
                       end
                   end
                   if ~intersection_possible
                       A_con = []; b_con = []; k_con = [];
                   end

                   obj.A_con_self{i}{j} = A_con;
                   obj.b_con_self{i}{j} = b_con;
                   obj.k_con_self{i}{j} = k_con;
               end
           end
        end
        
        function [h_obs, grad_h_obs] = evaluate_constraint(obj, A_con, b_con, k_con, k)
            n_k = size(k, 1);
            if (n_k ~= size(k_con, 1))
               error('number of rows of k does not match number of rows of k_con'); 
            end
            kappas = (k - obj.c_k(1:n_k))./obj.delta_k(1:n_k);
            
            % multiply rows of kappas together, replacing zeros with ones
            % (this yields the generator coefficients to be multiplied by
            % the fully-k-sliceable generators found in Vit_slc as in eqn.
            % (20) which have been combined into A_con)
            kappas_prod = double(k_con == 1).*kappas;
            kappas_prod(k_con ~= 1) = 1;
            kappas_prod = prod(kappas_prod, 1)';
            
            [h_obs_max, max_idx] = max(A_con*kappas_prod - b_con);
            h_obs = -h_obs_max;
            
            % compute gradients:
            % a little nasty, but the gradient will depend on the row of A
            % that produces the maximum h_obs_max, as well as which k's
            % that row depends on. if multiple rows of A create the same
            % maximum, then we take the maximum of the gradients (the
            % subgradient)
            if nargout == 2
                k_con_temp = (k_con == 1)';
                kappas_grad = double(k_con_temp);
                cols = 1:length(kappas);
                for i = cols
                    kappas_grad_temp = double(k_con_temp(:, i))*kappas(i);
                    kappas_grad_temp(k_con_temp(:, i) ~= 1) = 1;
                    kappas_grad(:, cols ~= i) = kappas_grad_temp.*kappas_grad(:, cols ~= i);
                end
                if length(max_idx) > 1
                    temp_grad_h_obs = A_con(max_idx, :)*kappas_grad;
                    grad_h_obs = (-max(temp_grad_h_obs)')./obj.delta_k(1:n_k);
                else
                    grad_h_obs = (-(A_con(max_idx, :)*kappas_grad)')./obj.delta_k(1:n_k);
                end
            end
        end
    end
end

