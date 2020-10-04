classdef three_link_rotatotope_FRS < robot_rotatotope_FRS_v2
    %THREE_LINK_ROTATOTOPE_FRS Represents the forward reachable set (FRS) of the
    %robot's links and joints as a set of rotatotopes.
    %   As described in Lemma 7, we overapproximate the FRS of a robot in
    %   workspace using sets of rotatotopes. We use these FRS's to
    %   construct collision-avoidance and self-intersection constraints,
    %   described in Sec. 4.B
    %
    
    properties

    end
    
    methods
        function obj = three_link_rotatotope_FRS(q_0, q_dot_0, FRS_options, varargin)
            
            % still to do: infer these properties from a particular agent
            % class, so that robot-specific subclasses are 
            % unnecessary:
            rotation_axes = [0 0 1 0 1 0;
                             0 1 0 1 0 1;
                             1 0 0 0 0 0]; % array of column vectors representing the body-fixed 3D rotation axis of each joint
            n_links = [3]; % number of physical robot links
            n_base = [0]; % number of "base links", not used for collision constraints but affect FRS construction
            n_joints = [6]; % number of joint for rotatotope_FRS
            link_predecessor_joints = {[1;2]; [1;2;3;4]; [1;2;3;4;5;6]}; % cell array of vectors, each specifying all predecessor joints of each link
            link_zonotopes = {zonotope([0.5, 0.5; 0, 0; 0, 0]);...
                              zonotope([0.5, 0.5; 0, 0; 0, 0]);...
                              zonotope([0.5, 0.5; 0, 0; 0, 0])}; % zonotopes representing the volume of each link
            joint_zonotopes = {zonotope([1; 0; 0]);...
                               zonotope([1; 0; 0]);...
                               zonotope([1; 0; 0])}; % zonotopes representing volume of joints at end of each link
            base_predecessor_joints = {}; % joints that precede base "link"... base link not used for collision detection, but affects FRS construction
            base_joint_zonotopes = {}; % represents "shoulder" base link of Fetch
            link_self_intersection = {[1;3]}; % cell array of links that could intersect for self-collision constraints
            
            
            pre_slice_dim = {[4]; [4]; [4]; [4]; [4]; [4]};  % for each JRS, dimensions to pre-slice
            pre_slice_values = {[q_dot_0(1)]; [q_dot_0(2)]; [q_dot_0(3)]; [q_dot_0(4)]; [q_dot_0(5)]; [q_dot_0(6)]}; % for each JRS, values to pre-slice

            k_dim = {[3]; [3]; [3]; [3]; [3]; [3]}; % for each JRS, dimensions of trajectory parameters
            k_names = {{'ka1'}; {'ka2'}; {'ka3'}; {'ka4'}; {'ka5'}; {'ka6'}}; % store list of all parameters
            
            % pass default args to super constructor
            obj@robot_rotatotope_FRS_v2(q_0, q_dot_0, FRS_options, ...
                'rotation_axes', rotation_axes, ...
                'n_links', n_links, ...
                'n_base', n_base, ...
                'n_joints', n_joints, ...
                'link_predecessor_joints', link_predecessor_joints, ...
                'link_zonotopes', link_zonotopes, ...
                'joint_zonotopes', joint_zonotopes, ...
                'base_predecessor_joints', base_predecessor_joints, ...
                'base_joint_zonotopes', base_joint_zonotopes, ...
                'link_self_intersection', link_self_intersection, ...
                'pre_slice_dim', pre_slice_dim, ...
                'pre_slice_values', pre_slice_values, ...
                'k_dim', k_dim, ...
                'k_names', k_names, ...
                varargin{:});
            
            % construct the FRS:
            obj = obj.create_FRS();
        end
    end
end

