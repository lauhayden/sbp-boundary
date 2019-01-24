function dydt = sbp_ipm(t, y, P, D1, u0_func, e0, u0_t_func, sigma)
    dydt = -P * D1 * (P * y + u0_func(t) * e0) + u0_t_func(t) * e0 - sigma * e0 * (y(1) - u0_func(t));
end