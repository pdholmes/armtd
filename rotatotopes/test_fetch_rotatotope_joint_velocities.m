clear all; clc;

figure(1); clf; hold on; axis equal;

%%% change this to try out pre-slicing!!
test_pre_slice = false;
pre_slice_interval = true;

% set FRS_options
FRS_options = struct();
FRS_options.buffer_dist = 0.001;
FRS_options.combs = generate_combinations_upto(200);
FRS_options.maxcombs = 200;
FRS_options.origin_shift = zeros(3, 1);

% q_0 = [0; 0; 0; pi/2 + 5*pi/48; 0; pi/2 + 5*pi/48];
q_0 = [-pi/2;0;0;0;0;0];
% q_dot_0 = [pi;pi/8;-pi/7;pi/9;pi/10;pi];
q_dot_0 = [pi; 0; 0; 0; 0; pi/2];
k = [0;0;0;0;0;0];

% R = robot_arm_FRS_rotatotope_fetch(q_0, q_dot_0, FRS_options);
% R = R.generate_constraints(O);

if test_pre_slice
    pre_slice_dim = {[5]; [5;4]; [5;4]; [5;4]; [5;4]; [5]};  % for each JRS, dimensions to pre-slice
    if pre_slice_interval
        % just try using a +- epsilon on acceleration for joints 2-5...
        % interval format is [lower_bound, upper_bound]
        epsilon = 0.1;
        pre_slice_values = {[q_dot_0(1)];...
            [q_dot_0(2), q_dot_0(2); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(3), q_dot_0(3); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(4), q_dot_0(4); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(5), q_dot_0(5); 0 - epsilon, 0 + epsilon];...
            [q_dot_0(6)]}; % for each JRS, intervals to pre-slice over
    else
        pre_slice_values = {[q_dot_0(1)]; [q_dot_0(2); 0]; [q_dot_0(3); 0]; [q_dot_0(4); 0]; [q_dot_0(5); 0]; [q_dot_0(6)]}; % for each JRS, values to pre-slice
    end
    k_dim = {[4]; []; []; []; []; [4]}; % for each JRS, dimensions of trajectory parameters
    k_names = {{'ka1'}; {}; {}; {}; {}; {'ka6'}}; % store list of all parameters 

    R2 = three_link_rotatotope_FRS(q_0, q_dot_0, FRS_options, ...
        'pre_slice_dim', pre_slice_dim, ...
        'pre_slice_values', pre_slice_values, ...
        'k_dim', k_dim, ...
        'k_names', k_names, ...
        'JRS_path', 'rotatotopes/joint_reachable_sets_v2/');
else
    % pre_slice_dim = {[5]; [5]; [5]; [5]; [5]; [5]};  % for each JRS, dimensions to pre-slice
    pre_slice_dim = {[6]; [6]; [6]; [6]; [6]; [6]};  % for each JRS, dimensions to pre-slice
    pre_slice_values = {q_dot_0(1); q_dot_0(2); q_dot_0(3); q_dot_0(4); q_dot_0(5); q_dot_0(6)};
    % k_dim = {[4]; [4]; [4]; [4]; [4]; [4]}; % for each JRS, dimensions of trajectory parameters
    k_dim = {[5]; [5]; [5]; [5]; [5]; [5]}; % for each JRS, dimensions of trajectory parameters
    k_names = {{'ka1'}; {'ka2'}; {'ka3'}; {'ka4'}; {'ka5'}; {'ka6'}}; % store list of all parameters 
    R2 = three_link_rotatotope_FRS(q_0, q_dot_0, FRS_options, ...
    'pre_slice_dim', pre_slice_dim, ...
    'pre_slice_values', pre_slice_values, ...
    'k_dim', k_dim, ...
    'k_names', k_names, ...
    'JRS_path', 'rotatotopes/joint_reachable_sets_v3/');
end


plot(R2);
figure(2); clf; hold on;
plot_EE_velocity(R2);
axis equal;
% plot_slice(R2, k);

% woo... so let's actually try to compute the velocity of the last 
% joint for some value of k, then see if that lies within the sliced rotatotope at that time

test_k = [-pi/3; pi/25; -pi/25; pi/25; pi/25; pi/6];
figure(5); clf; hold on;
% figure(2); hold on;
for i = 10:10:100
    R2.joint_velocity_rotatotopes{3}{i}.plot();
    R2.joint_velocity_rotatotopes{3}{i}.plot_slice(k_names, test_k);

    % plot joint positions:
    % R2.joint_position_rotatotopes{3}{i}.plot();
    % R2.joint_position_rotatotopes{3}{i}.plot_slice(k_names, test_k);
    %%% cool... verified the joint positions are correct. so something's just wrong with velocities.

    %...just gonna do this the cheap way, compute two positions and numerically differentiate:
    t1 = i/100 - 0.006;
    t2 = i/100 - 0.005;

    if t2 < 0.5
        q1 = q_0 +q_dot_0.*t1 + (1/2)*test_k.*t1.^2;
        q2 = q_0 +q_dot_0.*t2 + (1/2)*test_k.*t2.^2;
    else
        q_pk = q_0 +q_dot_0.*0.5 + (1/2)*test_k.*0.5.^2;
        q_dot_pk = q_dot_0 + test_k.*0.5;
        q1 = q_pk + q_dot_pk.*(t1 - 0.5) + 1/2*((0 - q_dot_pk)./0.5).*(t1 -0.5).^2;
        q2 = q_pk + q_dot_pk.*(t2 - 0.5) + 1/2*((0 - q_dot_pk)./0.5).*(t2 -0.5).^2;
    end

    % A = robot_arm_3D_fetch();
    A = arm_3D_3link;
    p1 = A.get_end_effector_location(q1);
    p2 = A.get_end_effector_location(q2);
    rotation_axes = [3, 2, 1, 2, 1, 2];
    for j = 1:6
        O1{j} = make_orientation(q1(j), rotation_axes(j));
        O2{j} = make_orientation(q2(j), rotation_axes(j));
    end

    p1_1 = O1{1}*O1{2}*[1;0;0] + O1{1}*O1{2}*O1{3}*O1{4}*[1;0;0] + O1{1}*O1{2}*O1{3}*O1{4}*O1{5}*O1{6}*[1;0;0];
    p2_1 = O2{1}*O2{2}*[1;0;0] + O2{1}*O2{2}*O2{3}*O2{4}*[1;0;0] + O2{1}*O2{2}*O2{3}*O2{4}*O2{5}*O2{6}*[1;0;0]; 

    estimated_velocity = (p2-p1)./(0.001);

    plot3(estimated_velocity(1), estimated_velocity(2), estimated_velocity(3), 'k.', 'MarkerSize', 30);
    % plot3(p2(1), p2(2), p2(3), 'k.', 'MarkerSize', 30);
end