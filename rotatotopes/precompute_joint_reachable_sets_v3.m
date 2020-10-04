function [] = precompute_joint_reachable_sets_v3()
%PRECOMPUTE_JOINT_REACHABLE_SETS (Offline) use CORA 2018 to compute the
%reachable sets of sines and cosines of an angle, parameterized by
%trajectory parameters K.
%   Section 3 (Offline Reachability Analysis) of the paper is implemented
%   here. First, we write down dynamics of cosine and sine of a joint angle
%   q_i that are parameterized by trajectory parameters K_i. Second, we
%   loop over subsets of the K^v_i dimension, computing reachable sets for
%   each. Finally, we save these zonotope joint reachable sets (JRSs) to be
%   used online by ARMTD.

% we're using the following state vector and parameterization:
% x = [cos(q_i); sin(q_i); k^a_i; k^v_i; t]
% where k^a_i is a desired acceleration on the interval t \in [0, t_plan]
% and k^v_i is the initial velocity of the joint.
%
% Including t (time) as a state isn't strictly necessary, but we will use
% it to precisely define when we switch to the failsafe maneuever of
% braking to a stop.
%
% UPDATED: 20201002: adding angular velocity as a state to keep track of
% within the JRSs.

% define hyperparameters for our particular parameterization:
dim = 7;
t_plan = 0.5;
t_total = 1;
dt = 0.01;

 % generates dynamics parameterized by K
[dyn_zero_to_t_plan, dyn_t_plan_to_t_total] = ...
    generate_parameterized_dynamics_v3(t_plan, t_total);

% described in Sec. 5.5: grid K^v_i space into smaller subintervals,
% compute separate JRSs in each one
n_JRS = 401; % separate initial velocity space (K^v_i) into 401 smaller intervals
c_kvi = linspace(-pi, pi, n_JRS); % centers of initial velocity subintervals
delta_kvi = (c_kvi(2) - c_kvi(1))/2; % subinterval is c_kvi +- delta_kvi

c_kai = 0; % acceleration parameter space (K^a_i) for each JRS centered at 0

% create folder to save precomputed JRSs
if ~exist('rotatotopes/joint_reachable_sets_v3', 'dir')
    mkdir('rotatotopes/joint_reachable_sets_v3');
end

% save vector of initial velocity subinterval centers
save('rotatotopes/joint_reachable_sets_v3/c_kvi.mat', 'c_kvi');

% set options for reachability analysis:
options.timeStep = dt;
options.taylorTerms=5; % number of taylor terms for reachable sets
options.zonotopeOrder= 2; % zonotope order... increase this for more complicated systems.
options.maxError = 1000*ones(dim, 1); % our zonotopes shouldn't be "splitting", so this term doesn't matter for now
options.verbose = 0;
options.uTrans = 0; % we won't be using any inputs, as traj. params specify trajectories
options.U = zonotope([0, 0]);
options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';

for j = 1:n_JRS % compute JRS for each velocity subinterval
    tic;
    % described in Sec. 5.5, scale K^a_i interval size with initial
    % velocity magnitude, lower bounded by pi/24 rad/s^2
    delta_kai = max(pi/24, abs(c_kvi(j)/3));
    
    % break JRS computation into two steps...
    % first, use dyn_zero_to_t_plan dynamics
    options.tStart = 0; % start time
    options.tFinal = t_plan; % end time for these dynamics
    
    % create initial zonotope for reachability analysis
    % this is described by eqs. (3) and (8)
    options.x0 = [1; 0; 0; c_kvi(j); c_kai; c_kvi(j); 0]; % start at q_i = 0, so cos(q_i) = 1 and sin(q_i) = 0
    % use two generators, one for K^a_i and one for K^v_i (eq. 8)
    options.R0 = zonotope([options.x0, [0; 0; 0; 0; delta_kai; 0; 0], [0; 0; 0; delta_kvi; 0; delta_kvi; 0]]);
    
    % create system for reachability analysis (1 dummy input)
    sys = nonlinearSys(dim, 1, dyn_zero_to_t_plan, options);
    
    % compute JRS over t \in [0, t_plan]
    JRS_zero_to_t_plan = reach(sys, options);
    
    % we'll use the last zonotope of JRS_zero_to_t_plan as the initial
    % zonotope for the braking dynamics portion of the JRS. however, we
    % will use the subset of the zonotope corresponding to t = t_plan:
    % (slicing described in Alg. 1)
    t_plan_slice = zonotope_slice(JRS_zero_to_t_plan{end}{1}, 7, t_plan);
    options.R0 = t_plan_slice;
    options.tStart = t_plan;
    options.tFinal = t_total;
    
    % create system for reachability analysis (1 dummy input)
    sys = nonlinearSys(dim, 1, dyn_t_plan_to_t_total, options);
    
    % compute JRS over t \in [t_plan, t_total]
    JRS_t_plan_to_t_total = reach(sys, options);
    
    % concatenate the full JRS (eq. 9, shown in Fig. 2)
    JRS = [JRS_zero_to_t_plan; JRS_t_plan_to_t_total];
    
    % save this JRS
    current_c_kvi = c_kvi(j);
    filename = sprintf('rotatotopes/joint_reachable_sets_v3/JRS_%0.3f.mat', current_c_kvi);
    save(filename, 'JRS', 'options', 't_plan', 't_total', 'current_c_kvi');
    
    % display:
    fprintf('Current initial velocity: %0.3f, ', current_c_kvi);
    toc;
end

