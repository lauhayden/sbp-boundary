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
sbp_operator = @D1_6;
solver = @rk4_wrapper;
%h = 0.0001;

% construct H for H-norm
H = cell(length(n), 1);
for i=1:length(n)
    [H_n, D1] = sbp_operator(n(i) + 1);
    H{i} = H_n / n(i);
end

% make cell arrays
% computed solutions
u_sat = cell(length(n), 1);
u_proj = cell(length(n), 1);
u_ipm = cell(length(n), 1);
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
    if igauss_on
        u_init_l = exp(-((x{i} - igauss_center) / igauss_sigma).^2);
    else
        u_init_l = zeros(size(x{i}));
    end
    u_init(i) = {u_init_l};
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
    [t_sat_l, y_sat_l] = linadv_solve(solver, 'sat', n_l, tf, u_init_l, ...
        sbp_operator, input_boundary_l, input_boundary_t_l, n_l); % TODO: CHECK SIGMA VAL
    [t_proj_l, y_proj_l] = linadv_solve(solver, 'proj', n_l, tf, u_init_l, ...
        sbp_operator, input_boundary_l, input_boundary_t_l, n_l);
    [t_ipm_l, y_ipm_l] = linadv_solve(solver, 'ipm', n_l, tf, u_init_l, ...
        sbp_operator, input_boundary_l, input_boundary_t_l, n_l);
    
    % store in cells
    u_sat(i) = {y_sat_l};
    u_proj(i) = {y_proj_l};
    u_ipm(i) = {y_ipm_l};
end

% calc errors
for i=1:length(n)
    % local vars
    solution_l = solution{i};
    
    % errors
    e_sat_l = solution_l - u_sat{i}(end, :)';
    e_proj_l = solution_l - u_proj{i}(end, :)';
    e_ipm_l = solution_l - u_ipm{i}(end, :)';
    
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
plot(x{1}, u_sat{1}(1, :)')
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
xlabel('x')
subplot(2, 2, 2)
plot(x{1}, u_sat{1}(end, :)')
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
xlabel('x')
subplot(2, 2, 3)
plot(x{1}, u_proj{1}(end, :)')
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
xlabel('x')
subplot(2, 2, 4)
plot(x{1}, u_ipm{1}(end, :)')
axis([0, 1, -0.1, 1.1])
title('SBP-IPM')
xlabel('x')

% plot error norms
figure
subplot(2, 2, 2)
loglog(n, enorm_sat, 'o-')
title('SBP-SAT Error')
ylabel('H-norm of error')
xlabel('1/h')
subplot(2, 2, 3)
loglog(n, enorm_proj, 'o-')
title('SBP-Proj Error')
ylabel('H-norm of error')
xlabel('1/h')
subplot(2, 2, 4)
loglog(n, enorm_ipm, 'o-')
title('SBP-IPM Error')
ylabel('H-norm of error')
xlabel('1/h')

% plot exact solution overlaid
% figure
% hold on
% plot(x{4}, solution{4}, 'o')
% plot(x{4}, u_sat{4}(end, :)', 'o')
% legend('exact', 'sat4')

% plot errors
figure
hold on
for i=1:length(n)
    plot(x{i}, solution{i} - u_sat{i}(end, :)')
end
title('Error at tf')
ylabel('u\_exact - u\_sbp\_sat')
xlabel('x')

figure
hold on
for i=1:length(n)
    plot(x{i}, solution{i} - u_proj{i}(end, :)')
end
title('Error at tf')
ylabel('u\_exact - u\_sbp\_proj')
xlabel('x')

% plot boundary overlaid with exact
% t_boundary = linspace(0, 0.5, 501);
% y_boundary = input_boundary(t_boundary, true, bgauss_center, bgauss_sigma);
% figure
% hold on
% plot(x{4}(1:501), solution{4}(1:501) - y_boundary')

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