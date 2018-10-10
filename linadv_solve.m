function [t, u] = linadv_solve(n, tf, u_init, u0_func, D1_func)
    [H, D1] = D1_func(n);
    
    % scale SBP operator to grid size
    H = H / n;
    D1 = D1 * n;
    
    % precalculate element 1, 1 of inv(H)
    invH = inv(H);
    invH_11 = invH(1, 1);

    function dydt = sbp_sat(t, y, n, u0_func, invH_11, D1)
        e0 = zeros(n, 1);
        e0(1) = 1;
        dydt = -(D1 * y) - (y(1) - u0_func(t)) * invH_11 * e0;
    end

    [t, u] = ode45(@(t, y) sbp_sat(t, y, n, u0_func, invH_11, D1), [0, tf], u_init);
end