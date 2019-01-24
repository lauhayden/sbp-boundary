function dydt = sbp_ipm(t, y, D1, u0_func, u0_t_func, sigma)
    dydt = y;
    dydt(1) = u0_func(t);
    dydt = -D1 * dydt;
    dydt(1) = u0_t_func(t) - sigma * (y(1) - u0_func(t));
end