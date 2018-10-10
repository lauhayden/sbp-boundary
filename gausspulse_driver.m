n = 200;
tf = 0.5;
sigma = 0.1;
center = 0.5;

x = linspace(0, 1, n);
u_init = exp(-((x - center) / sigma).^2);

[t, y] = linadv_solve(n, tf, u_init, @(t) 0, @D1_6);

plot(x, y(end, :))
axis([0, 1, -0.1, 1.1])