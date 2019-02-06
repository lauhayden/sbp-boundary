% inconsistent boundary data driver
% Initial state is all 1, homogenous BCs

close all

n = 200;
tf = 0.5;
h = 1 / n;


x = linspace(0, 1, n + 1);
u_init = ones(n + 1, 1);

[t, y1] = linadv_solve(@ode45, 'sat', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);
[t, y2] = linadv_solve(@ode45, 'proj', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);
[t, y3] = linadv_solve(@ode45, 'ipm', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);

subplot(2, 2, 1)
plot(x, y1(1, :))
axis([0, 1, -0.1, 1.1])
title('Initial Condition')
xlabel('x')
subplot(2, 2, 2)
plot(x, y1(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-SAT')
xlabel('x')
subplot(2, 2, 3)
plot(x, y2(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-Proj')
xlabel('x')
subplot(2, 2, 4)
plot(x, y3(end, :))
axis([0, 1, -0.1, 1.1])
title('SBP-IPM')
xlabel('x')