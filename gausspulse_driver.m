% Simple double gaussian pulse driver for linear advection equation
% Initial state includes one gaussian, boundary condition adds another

close all

% setup params
n = [100, 200, 500, 1000];
tf = 0.5;
igauss_on = true;
igauss_sigma = 0.05;
igauss_center = 0.25;
bgauss_on = true;
bgauss_sigma = 0.05;
bgauss_center = -0.25;

% construct H for H-norm
H = cell(length(n), 1);
for i=1:length(n)
    [H_n, D1] = D1_6(n(i) + 1);
    H{i} = H_n / n(i);
end

% make cell arrays
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
% error norms
enorm_sat = zeros(length(n), 1);
enorm_proj = zeros(length(n), 1);
enorm_ipm = zeros(length(n), 1);

% x-values
x = cell(length(n), 1);
for i=1:length(n)
    x(i) = {linspace(0, 1, n(i) + 1)'};
end

% initial condition
u_init = cell(length(n), 1);
for i=1:length(n)
    u_init(i) = {exp(-((x{i} - igauss_center) / igauss_sigma).^2)};
end

% exact solution
solution = cell(length(n), 1);
for i=1:length(n)
    solution_l = zeros(n(i) + 1, 1);
    if igauss_on
        solution_l = solution_l + exp(-((x{i} - igauss_center - tf) / igauss_sigma).^2);
    end
    if bgauss_on
        solution_l = solution_l + exp(-((x{i} - bgauss_center - tf) / bgauss_sigma).^2);
    end
    solution(i) = {solution_l};
end
    
%% Computation

for i=1:length(n)
    % local vars
    n_l = n(i);
    u_init_l = u_init{i};
    input_boundary_l = @(t) input_boundary(t, bgauss_on, bgauss_center, bgauss_sigma);
    input_boundary_t_l = @(t) input_boundary_t(t, bgauss_on, bgauss_center, bgauss_sigma);
    
    % solve
    [t_sat_l, y_sat_l] = linadv_solve('sat', n_l, tf, u_init_l, ...
        @D1_6, input_boundary_l, input_boundary_t_l, n_l); % CHECK SIGMA VAL
    [t_proj_l, y_proj_l] = linadv_solve('proj', n_l, tf, u_init_l, ...
        @D1_6, input_boundary_l, input_boundary_t_l, n_l);
    [t_ipm_l, y_ipm_l] = linadv_solve('ipm', n_l, tf, u_init_l, ...
        @D1_6, input_boundary_l, input_boundary_t_l, n_l);
    
    % store in cells
    y_sat(i) = {y_sat_l};
    y_proj(i) = {y_proj_l};
    y_ipm(i) = {y_ipm_l};
end

% calc errors
for i=1:length(n)
    % local vars
    solution_l = solution{i};
    
    % errors
    e_sat_l = solution_l - y_sat{i}(end, :)';
    e_proj_l = solution_l - y_proj{i}(end, :)';
    e_ipm_l = solution_l - y_ipm{i}(end, :)';
    
    % error norms
    enorm_sat(i) = sqrt(e_sat_l' * H{i} * e_sat_l);
    enorm_proj(i) = sqrt(e_proj_l' * H{i} * e_proj_l);
    enorm_ipm(i) = sqrt(e_ipm_l' * H{i} * e_ipm_l);
    
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
plot(x{1}, y_sat{1}(1, :)')
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
subplot(2, 2, 2)
plot(x{1}, y_sat{1}(end, :)')
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
subplot(2, 2, 3)
plot(x{1}, y_proj{1}(end, :)')
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
subplot(2, 2, 4)
plot(x{1}, y_ipm{1}(end, :)')
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
plot(x{4}, solution{4}, 'o')
plot(x{4}, y_sat{4}(end, :)', 'o')

% plot errors
figure
hold on
for i=1:length(n)
    plot(x{i}, solution{i} - y_sat{i}(end, :)')
end

%% Inline functions

function x0 = input_boundary(t, bgauss_on, bgauss_center, bgauss_sigma)
    if bgauss_on
        x0 = exp(-((t + bgauss_center) / bgauss_sigma).^2);
    else
        x0 = 0;
    end
end

function x0_t = input_boundary_t(t, bgauss_on, bgauss_center, bgauss_sigma)
    if bgauss_on
        x0_t = exp(-((t + bgauss_center) / bgauss_sigma).^2);
        x0_t = x0_t * -2 / bgauss_sigma^2 * (t + bgauss_center);
    else
        x0_t = 0;
    end
end