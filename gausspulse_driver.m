% Simple double gaussian pulse driver for linear advection equation
% Initial state includes one gaussian, boundary condition adds another

close all

n_values = [100, 200, 500, 1000];
tf = 0.5;
sigma = 0.05;
center = 0.25;

H = cell(length(n_values), 1);

for i=1:length(n_values)
    [H_n, D1] = D1_6(n_values(i));
    H{i} = H_n / n_values(i);
end

y_sat = cell(length(n_values), 1);
y_proj = cell(length(n_values), 1);
y_ipm = cell(length(n_values), 1);

for i=1:length(n_values)
    n = n_values(i);
    x = linspace(0, 1, n);
    u_init = exp(-((x - center) / sigma).^2);

    [t, y1] = linadv_solve('sat', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n);
    [t, y2] = linadv_solve('proj', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n);
    [t, y3] = linadv_solve('ipm', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, n);
    y_sat(i) = {y1};
    y_proj(i) = {y2};
    y_ipm(i) = {y3};
end

subplot(2, 2, 1)
plot(linspace(0, 1, n_values(1)), y_sat{1}(1, :))
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
subplot(2, 2, 2)
plot(linspace(0, 1, n_values(1)), y_sat{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
subplot(2, 2, 3)
plot(linspace(0, 1, n_values(1)), y_proj{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
subplot(2, 2, 4)
plot(linspace(0, 1, n_values(1)), y_ipm{1}(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-IPM')

error_norms_sat = zeros(length(n_values), 1);
error_norms_proj = zeros(length(n_values), 1);
error_norms_ipm = zeros(length(n_values), 1);

for i=1:length(n_values)
    n = n_values(i);
    solution = exact_solution(tf, n);
    error_sat = solution - y_sat{i}(end, :);
    error_proj = solution - y_proj{i}(end, :);
    error_ipm = solution - y_ipm{i}(end, :);
    error_norms_sat(i) = sqrt(error_sat * H{i} * error_sat');
    error_norms_proj(i) = sqrt(error_proj * H{i} * error_proj');
    error_norms_ipm(i) = sqrt(error_ipm * H{i} * error_ipm');
end

figure
subplot(2, 2, 2)
loglog(n_values, error_norms_sat, 'o-')
title('SBP-SAT Error')
subplot(2, 2, 3)
loglog(n_values, error_norms_proj, 'o-')
title('SBP-Proj Error')
subplot(2, 2, 4)
loglog(n_values, error_norms_ipm, 'o-')
title('SBP-IPM Error')

A_slope = ones(length(n_values), 2);
A_slope(:, 2) = log10(n_values);
regression_sat = A_slope \ log10(error_norms_sat);
regression_proj = A_slope \ log10(error_norms_proj);
regression_ipm = A_slope \ log10(error_norms_ipm);

disp(regression_sat)
disp(regression_proj)
disp(regression_ipm)

figure
hold on
plot(linspace(0, 1, 200), exact_solution(tf, 200))
plot(linspace(0, 1, 200), y_sat{2}(end, :))

function soln = exact_solution(t, n)
    x = linspace(0, 1, n);
    soln = exp(-((x - 0.25 - t) / 0.05).^2);
    soln = soln + exp(-((x + 0.25 - t) / 0.05).^2);
end

function x0 = input_boundary(t)
    x0_center = 0.25;
    x0_sigma = 0.05;
    x0 = exp(-((t - x0_center) / x0_sigma).^2);
end

function x0_t = input_boundary_t(t)
    x0_center = 0.25;
    x0_sigma = 0.05;
    x0_t = exp(-((t - x0_center) / x0_sigma).^2);
    x0_t = x0_t * -2 / x0_sigma^2 * (t - x0_center);    
end