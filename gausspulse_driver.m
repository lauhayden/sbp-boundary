% Simple double gaussian pulse driver for linear advection equation
% Initial state includes one gaussian, boundary condition adds another

close all

% setup params
n = [100, 200, 500, 1000];
tf = 0.25;
sigma = 0.05;
center = 0.5;

% construct H for H-norm
H = cell(length(n), 1);
for i=1:length(n)
    [H_n, D1] = D1_6(n(i) + 1);
    H{i} = H_n / n(i);
end

% make cell arrays
% x-values
x = cell(length(n), 1);
% computed solutions
y_sat = cell(length(n), 1);
y_proj = cell(length(n), 1);
y_ipm = cell(length(n), 1);
% computed times
t_sat = cell(length(n), 1);
t_proj = cell(length(n), 1);
t_ipm = cell(length(n), 1);
% errors at end time
e_sat = cell(length(n), 1);
e_proj = cell(length(n), 1);
e_ipm = cell(length(n), 1);
enorm_sat = zeros(length(n), 1);
enorm_proj = zeros(length(n), 1);
enorm_ipm = zeros(length(n), 1);

%% Computation

for i=1:length(n)
    % local vars
    n_l = n(i);
    x_l = linspace(0, 1, n_l + 1);
    u_init = exp(-((x_l - center) / sigma).^2);

    % solve
    [t_sat_l, y_sat_l] = linadv_solve('sat', n_l, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n_l); % CHECK SIGMA VAL
    [t_proj_l, y_proj_l] = linadv_solve('proj', n_l, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n_l);
    [t_ipm_l, y_ipm_l] = linadv_solve('ipm', n_l, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n_l);
    
    % store in cells
    x(i) = {x_l};
    y_sat(i) = {y_sat_l};
    y_proj(i) = {y_proj_l};
    y_ipm(i) = {y_ipm_l};
end

% calc errors
for i=1:length(n)
    % local vars
    x_l = n(i);
    solution = exact_solution(tf, x_l);
    
    % errors
    e_sat_l = solution - y_sat{i}(end, :);
    e_proj_l = solution - y_proj{i}(end, :);
    e_ipm_l = solution - y_ipm{i}(end, :);
    
    % error norms
    enorm_sat(i) = sqrt(e_sat_l * H{i} * e_sat_l');
    enorm_proj(i) = sqrt(e_proj_l * H{i} * e_proj_l');
    enorm_ipm(i) = sqrt(e_ipm_l * H{i} * e_ipm_l');
    
    % store in cells
    e_sat(i) = {e_sat_l};
    e_proj(i) = {e_proj_l};
    e_ipm(i) = {e_ipm_l};
end

% linear regression on log-log plot of error norms
A_slope = ones(length(n), 2);
A_slope(:, 2) = log10(n);
enormreg_sat = A_slope \ log10(enorm_sat);
enormreg_proj = A_slope \ log10(enorm_proj);
enormreg_ipm = A_slope \ log10(enorm_ipm);

%% Plotting

% plot final solutions together
subplot(2, 2, 1)
plot(x{1}, y_sat{1}(1, :))
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
subplot(2, 2, 2)
plot(x{1}, y_sat{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
subplot(2, 2, 3)
plot(x{1}, y_proj{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
subplot(2, 2, 4)
plot(x{1}, y_ipm{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-IPM')

% plot error norms
figure
subplot(2, 2, 2)
loglog(n, enorm_sat, 'o-')
title('SBP-SAT Error')
subplot(2, 2, 3)
loglog(n, enorm_proj, 'o-')
title('SBP-Proj Error')
subplot(2, 2, 4)
loglog(n, enorm_ipm, 'o-')
title('SBP-IPM Error')

% plot exact solution overlaid
figure
hold on
plot(x{4}, exact_solution(tf, x{4}), 'o')
plot(x{4}, y_sat{4}(end, :), 'o')

% plot errors
figure
hold on
for i=1:length(n)
    plot(x{i}, exact_solution(tf, x{i}) - y_sat{i}(end, :))
end

% figure
% hold on
% plot(linspace(0, 1, n_values(1)), exact_solution(0, n_values(1)))
% plot(linspace(0, 1, n_values(1)), y_sat{1}(1, :))

%% Inline functions

function soln = exact_solution(t, x)
    soln = exp(-((x - 0.5 - t) / 0.05).^2);
%     soln = soln + exp(-((x + 0.25 - t) / 0.05).^2);
end

function x0 = input_boundary(t)
%     x0_center = 0.25;
%     x0_sigma = 0.05;
%     x0 = exp(-((t - x0_center) / x0_sigma).^2);
    x0 = 0;
end

function x0_t = input_boundary_t(t)
%     x0_center = 0.25;
%     x0_sigma = 0.05;
%     x0_t = exp(-((t - x0_center) / x0_sigma).^2);
%     x0_t = x0_t * -2 / x0_sigma^2 * (t - x0_center);    
    x0_t = 0;
end