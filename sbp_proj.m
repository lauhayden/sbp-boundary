function dydt = sbp_proj(t, y, D1, u0_t_func)
    dydt = -D1 * y;
    dydt(1) = u0_t_func(t);
end