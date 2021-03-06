clear all; clc;

figure(1); clf; hold on; axis equal;

% set FRS_options
FRS_options = struct();
FRS_options.buffer_dist = 0.1460;
FRS_options.combs = generate_combinations_upto(200);
FRS_options.maxcombs = 200;
FRS_options.origin_shift = zeros(3, 1);

obs_center = [0.3877; 1.018; -0.008];
obs_width = [0.1];
O{1} = box_obstacle_zonotope('center', obs_center(:), 'side_lengths', [obs_width, obs_width, obs_width]);
plot(O{1});

% q_0 = [0; 0; 0; pi/2 + 5*pi/48; 0; pi/2 + 5*pi/48];
q_0 = [-pi/2;0;0;0;0;0];
q_dot_0 = [pi;0;0;0;0;pi/2];
k = [0;0;0;0;0;0];

R = robot_arm_FRS_rotatotope_fetch(q_0, q_dot_0, FRS_options);
R = R.generate_constraints(O);
R2 = fetch_rotatotope_FRS(q_0, q_dot_0, FRS_options);
R2 = R2.generate_constraints_from_obstacles(O);

% you can check that R and R2's A_con, b_con, and k_con properties are
% identical!

plot(R2);
% plot_slice(R2, k);

% test the constraints at A_con{1}{3}{100} (final time step of 3rd link for
% 1st obstacle). we'll set k2 = k3 = k4 = k5 = 0, and grid over the k1, k6
% dimensions.

idx1 = 1;
idx2 = 3;
idx3 = 100;

A_con = R2.A_con{idx1}{idx2}{idx3};
b_con = R2.b_con{idx1}{idx2}{idx3};
k_con = R2.k_con{idx1}{idx2}{idx3};


c_k = R2.c_k;
delta_k = R2.delta_k;

figure(2); clf; hold on;
lims = [-delta_k(1) -delta_k(1) delta_k(1) delta_k(1) -delta_k(1); -delta_k(6) delta_k(6) delta_k(6) -delta_k(6) -delta_k(6)];
plot(lims(1, :)', lims(2, :)', 'k--', 'LineWidth', 4);

myk1 = linspace(-delta_k(1), delta_k(1), 100);
myk6 = linspace(-delta_k(6), delta_k(6), 100);
[K1_grid, K6_grid] = meshgrid(myk1, myk6);
h_obs = zeros(length(myk1),length(myk6));
for i = 1:length(myk1)
    for j = 1:length(myk6)
        k = [K1_grid(i, j); 0; 0; 0; 0; K6_grid(i, j)];
        h_obs(i, j) = R2.evaluate_constraint(A_con, b_con, k_con, k);
        if h_obs(i, j) >= 0
            plot(K1_grid(i, j), K6_grid(i, j), 'r.', 'MarkerSize', 10);
        end
    end
end

xlabel('K_1', 'FontSize', 20);
ylabel('K_6', 'FontSize', 20);
title('Red trajectory parameters could cause collision');

disp('Click a point!');
[p1, p2] = ginput(1);
plot(p1, p2, 'kx', 'MarkerSize', 16, 'LineWidth', 6);

figure(3); clf; axis equal; hold on;
slice_pt = [p1;0;0;0;0;p2];
R2.plot_slice(slice_pt, 100);
O{1}.plot_data.body = [];
plot(O{1});
