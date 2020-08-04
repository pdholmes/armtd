function [dyn_zero_to_t_plan, dyn_t_plan_to_t_total] = generate_parameterized_dynamics(t_plan,t_total)
%GENERATE_PARAMETERIZED_DYNAMICS prepares dynamics for CORA 2018 with the
%particular parameterization described in Section 2 of the paper.
%   [dyn_zero_to_t_plan, dyn_t_plan_to_t_total] = generate_parameterized_dynamics(t_plan, t_total)
%   starting from initial velocity k_i^v, we use constant acceleration over
%   t \in [0, t_plan]. 
%   then, we define trajectories with a failsafe
%   (braking) maneuver from the peak speed over t \in [t_plan, t_total].

if ~exist('rotatotopes/dynamics', 'dir')
   mkdir('rotatotopes/dynamics'); 
end

syms cqi sqi kai kvi t real;
syms tdummy udummy real; % CORA will require these arguments, but we won't use them.
x = [cqi; sqi; kai; kvi; t];

% these dynamics are written in eqs. (2) and (5) in the paper
q_i_dot = kvi + kai*t;
dcqi = -sqi*q_i_dot;
dsqi = cqi*q_i_dot;
dkai = 0;
dkvi = 0;
dt = 1;

dx = [dcqi; dsqi; dkai; dkvi; dt];
dyn_zero_to_t_plan = matlabFunction(dx, 'File', 'rotatotopes/dynamics/dyn_zero_to_t_plan', 'vars', {tdummy, x, udummy});

% now we specify braking dynamics on t \in [t_plan, t_total]
t_to_stop = t_total - t_plan;
q_i_dot_pk = kvi + kai*t_plan;
braking_acceleration = (0 - q_i_dot_pk)/t_to_stop; % brake to 0 velocity from q_i_dot_pk in t_to_stop seconds
q_i_dot = q_i_dot_pk + braking_acceleration*(t - t_plan);

dcqi = -sqi*q_i_dot;
dsqi = cqi*q_i_dot;

dx = [dcqi; dsqi; dkai; dkvi; dt];
dyn_t_plan_to_t_total = matlabFunction(dx, 'File', 'rotatotopes/dynamics/dyn_t_plan_to_t_total', 'vars', {tdummy, x, udummy});

end

