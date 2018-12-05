% inconsistent boundary data driver
% Initial state is all 1, homogenous BCs

close all

n = 200;
tf = 0.5;
h = 1 / n;

x = linspace(0, 1, n);
u_init = ones(1, n);

[t, y1] = linadv_solve('sat', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);
[t, y2] = linadv_solve('proj', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);
[t, y3] = linadv_solve('ipm', n, tf, u_init, @D1_6, @(t) 0, @(t) 0, 1/h);

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