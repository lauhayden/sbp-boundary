function dydt = sbp_sat(t, y, D1, u0_func, invH_11, e0)
    dydt = -(D1 * y) - (y(1) - u0_func(t)) * invH_11 * e0;
end