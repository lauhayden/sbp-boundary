function dydt = sbp_sat(t, y, D1, u0_func, invH_11)
    dydt = -D1 * y;
    dydt(1) = dydt(1) - (y(1) - u0_func(t)) * invH_11;
end