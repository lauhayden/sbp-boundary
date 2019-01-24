function dydt = sbp_proj(t, y, P, D1, u0_t_func, e0)
    dydt = -P * D1 * y + u0_t_func(t) * e0;
end