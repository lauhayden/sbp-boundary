% Simple double gaussian pulse driver for linear advection equation
% Initial state includes one gaussian, boundary condition adds another

close all

n = 200;
tf = 0.5;
sigma = 0.1;
center = 0.3;
h = 1 / n;

x = linspace(0, 1, n);
u_init = exp(-((x - center) / sigma).^2);

[t, y1] = linadv_solve('sat', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, 1/h);
[t, y2] = linadv_solve('proj', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, 1/h);
[t, y3] = linadv_solve('ipm', n, tf, u_init, @D1_6, @input_boundary, @input_boundary_t, 1/h);

subplot(2, 2, 1)
plot(x, y1(1, :))
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
subplot(2, 2, 2)
plot(x, y1(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
subplot(2, 2, 3)
plot(x, y2(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
subplot(2, 2, 4)
plot(x, y3(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-IPM')

function x0 = input_boundary(t)
    x0_center = 0.25;
    x0_sigma = 0.05;
    if t < x0_center - 3 * x0_sigma || t > x0_center + 3 * x0_sigma
        x0 = 0;
    else
        x0 = exp(-((t - x0_center) / x0_sigma).^2);
    end
end

function x0_t = input_boundary_t(t)
    x0_center = 0.25;
    x0_sigma = 0.05;
    if t < x0_center - 3 * x0_sigma || t > x0_center + 3 * x0_sigma
        x0_t = 0;
    else
        x0_t = exp(-((t - x0_center) / x0_sigma).^2);
        x0_t = x0_t * -2 / x0_sigma^2 * (t - x0_center);    
    end
end