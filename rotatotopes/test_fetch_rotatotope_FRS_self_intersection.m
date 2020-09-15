clear all; clc;

figure(1); clf; hold on; axis equal;

%%% change this to try out pre-slicing!!
test_pre_slice = true;
pre_slice_interval = true;

% set FRS_options
FRS_options = struct();
FRS_options.buffer_dist = 0.1460;
FRS_options.combs = generate_combinations_upto(200);
FRS_options.maxcombs = 200;
FRS_options.origin_shift = zeros(3, 1);

q_0 = [0; 0; 0; pi/4; 0; pi/4];
q_dot_0 = [0; 0; 0; pi/2; 0; pi/4 + 0.1];

% R = robot_arm_FRS_rotatotope_fetch(q_0, q_dot_0, FRS_options);
% R = R.generate_self_intersection_constraints();

if test_pre_slice
    pre_slice_dim = {[4;3]; [4;3]; [4;3]; [4]; [4;3]; [4]};  % for each JRS, dimensions to pre-slice
    if pre_slice_interval
        % just try using a +- epsilon on acceleration for joints 2-5...
        % interval format is [lower_bound, upper_bound]
        epsilon = 0.1;
        pre_slice_values = {[q_dot_0(1), q_dot_0(1); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(2), q_dot_0(2); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(3), q_dot_0(3); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(4)];...
            [q_dot_0(5), q_dot_0(5); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(6)]}; % for each JRS, intervals to pre-slice over
    else
        pre_slice_values = {[q_dot_0(1); 0]; [q_dot_0(2); 0]; [q_dot_0(3); 0]; [q_dot_0(4)]; [q_dot_0(5); 0]; [q_dot_0(6)]}; % for each JRS, values to pre-slice
    end
    k_dim = {[]; []; []; [3]; []; [3]}; % for each JRS, dimensions of trajectory parameters
    k_names = {{}; {}; {}; {'ka4'}; {}; {'ka6'}}; % store list of all parameters 

    R2 = fetch_rotatotope_FRS(q_0, q_dot_0, FRS_options, ...
        'pre_slice_dim', pre_slice_dim, ...
        'pre_slice_values', pre_slice_values, ...
        'k_dim', k_dim, ...
        'k_names', k_names);
else
    R2 = fetch_rotatotope_FRS(q_0, q_dot_0, FRS_options);
end

R2 = R2.generate_self_intersection_constraints();

% you can check that R and R2's A_con, b_con, and k_con properties are
% identical! (in the test_pre_slice = false case)

plot(R2);
% plot_slice(R2, k);

% test the self-intersection constraints at A_con_self{1}{100} (final time
% step for possible 1st and 3rd link intersection). 
% we'll set k1 = k2 = k3 = k5 = 0, and grid over the k4, k6 dimensions.

idx1 = 1;
idx2 = 100;

A_con_self = R2.A_con_self{idx1}{idx2};
b_con_self = R2.b_con_self{idx1}{idx2};
k_con_self = R2.k_con_self{idx1}{idx2};

c_k = R2.c_k;
delta_k = R2.delta_k;

figure(2); clf; hold on;
if test_pre_slice
    lims = [-delta_k(1) -delta_k(1) delta_k(1) delta_k(1) -delta_k(1); -delta_k(2) delta_k(2) delta_k(2) -delta_k(2) -delta_k(2)];
    plot(lims(1, :)', lims(2, :)', 'k--', 'LineWidth', 4);
    myk4 = linspace(-delta_k(1), delta_k(1), 100);
    myk6 = linspace(-delta_k(2), delta_k(2), 100);
else
    lims = [-delta_k(4) -delta_k(4) delta_k(4) delta_k(4) -delta_k(4); -delta_k(6) delta_k(6) delta_k(6) -delta_k(6) -delta_k(6)];
    plot(lims(1, :)', lims(2, :)', 'k--', 'LineWidth', 4);
    myk4 = linspace(-delta_k(4), delta_k(4), 100);
    myk6 = linspace(-delta_k(6), delta_k(6), 100);
end

[K4_grid, K6_grid] = meshgrid(myk4, myk6);
h_obs = zeros(length(myk4),length(myk6));
for i = 1:length(myk4)
    for j = 1:length(myk6)
        if test_pre_slice
            k = [K4_grid(i, j); K6_grid(i, j)];
            h_obs(i, j) = R2.evaluate_constraint(A_con_self, b_con_self, k_con_self, {'ka4'; 'ka6'}, k);
        else
            k = [0; 0; 0; K4_grid(i, j); 0; K6_grid(i, j)];
            h_obs(i, j) = R2.evaluate_constraint(A_con_self, b_con_self, k_con_self, {'ka1'; 'ka2'; 'ka3'; 'ka4'; 'ka5'; 'ka6'}, k);
        end
        if h_obs(i, j) >= 0
            plot(K4_grid(i, j), K6_grid(i, j), 'r.', 'MarkerSize', 10);
        end
    end
end

xlabel('K_4', 'FontSize', 20);
ylabel('K_6', 'FontSize', 20);
title('Red trajectory parameters could cause self-intersection');

disp('Click a point!');
[p1, p2] = ginput(1);
plot(p1, p2, 'kx', 'MarkerSize', 16, 'LineWidth', 6);

figure(3); clf; axis equal; hold on;
if test_pre_slice
    slice_pt = [p1;p2];
    R2.plot_slice({'ka4'; 'ka6'}, slice_pt, 100);
else
    slice_pt = [0;0;0;p1;0;p2];
    R2.plot_slice({'ka1'; 'ka2'; 'ka3'; 'ka4'; 'ka5'; 'ka6'}, slice_pt, 100);
end