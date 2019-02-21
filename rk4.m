% RK4 integrator
function u_next = rk4(h, func, t, u)
    k1 = h * func(t, u);
    k2 = h * func(t + h / 2, u + k1 / 2);
    k3 = h * func(t + h / 2, u + k2 / 2);
    k4 = h * func(t + h, u + k3);
    u_next = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end