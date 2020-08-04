classdef fetch_rotatotope_FRS < robot_rotatotope_FRS
    %FETCH_ROTATOTOPE_FRS Represents the forward reachable set (FRS) of the
    %robot's links and joints as a set of rotatotopes.
    %   As described in Lemma 7, we overapproximate the FRS of a robot in
    %   workspace using sets of rotatotopes. We use these FRS's to
    %   construct collision-avoidance and self-intersection constraints,
    %   described in Sec. 4.B
    %
    %   This class is a subclass for Fetch containing
    %   robot-specific parameters.
    
    properties

    end
    
    methods
        function obj = fetch_rotatotope_FRS(varargin)
            obj = obj@robot_rotatotope_FRS(varargin{:});
            
            % still to do: infer these properties from a particular agent
            % class, so that robot-specific subclasses are 
            % unnecessary:
            obj.rotation_axes = [0 0 1 0 1 0;
                             0 1 0 1 0 1;
                             1 0 0 0 0 0]; % array of column vectors representing the body-fixed 3D rotation axis of each joint
            obj.n_links = [3]; % number of physical robot links
            obj.n_base = [1]; % number of "base links", not used for collision constraints but affect FRS construction
            obj.link_predecessor_joints = {[1;2], [1;2;3;4], [1;2;3;4;5;6]}; % cell array of vectors, each specifying all predecessor joints of each link
            obj.link_zonotopes = {zonotope([0.3556/2, 0.3556/2; 0, 0; 0, 0]);...
                              zonotope([0.3302/2, 0.3302/2; 0, 0; 0, 0]);...
                              zonotope([0.3302/2, 0.3302/2; 0, 0; 0, 0])}; % zonotopes representing the volume of each link
            obj.joint_zonotopes = {zonotope([0.3556; 0; 0]);...
                               zonotope([0.3302; 0; 0]);...
                               zonotope([0.3302; 0; 0])}; % zonotopes representing volume of joints at end of each link
            obj.base_predecessor_joints = {[1]}; % joints that precede base "link"... base link not used for collision detection, but affects FRS construction
            obj.base_joint_zonotopes = {zonotope([0.1206; 0; 0.0825])}; % represents "shoulder" base link of Fetch
            obj.link_self_intersection = {[1;3]}; % cell array of links that could intersect for self-collision constraints
            obj.k_dim = [3];
            
            % construct the FRS:
            obj = obj.create_FRS();
        end
    end
end

