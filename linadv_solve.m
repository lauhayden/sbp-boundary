% Simple Linear Advection Equation solver using SBP

% n : grid size [0:1]
% tf : final time
% u_init : initial condition
% D1_func : SBP operator function
% u0_func : boundary condition at x=0, function of time
% u0_t_func : time derivative of boundary condition at x=0

function [t, u] = linadv_solve(bc_method, n, tf, u_init, D1_func, u0_func, u0_t_func, sigma)
    [H, D1] = D1_func(n);
    
    % scale SBP operator to grid size
    H = H / n;
    D1 = D1 * n;
    
    % precalculate element 1, 1 of inv(H)
    invH = inv(H);
    invH_11 = invH(1, 1);
    
    % precalc other vectors and matrices
    e0 = zeros(n, 1);
    e0(1) = 1;
    L = e0;
    P = eye(n) - invH * L * inv(L' * invH * L) * L';
    
    function dydt = sbp_sat(t, y)
        dydt = -(D1 * y) - (y(1) - u0_func(t)) * invH_11 * e0;
    end

    function dydt = sbp_proj(t, y)
        dydt = -P * D1 * y + u0_t_func(t) * e0;
    end

    function dydt = sbp_ipm(t, y)
        dydt = -P * D1 * (P * y + u0_func(t) * e0) + u0_t_func(t) * e0 - sigma * e0 * (y(1) - u0_func(t));
    end

    if bc_method == "sat"
        [t, u] = ode45(@sbp_sat, [0, tf], u_init);
    elseif bc_method == "proj"
        [t, u] = ode45(@sbp_proj, [0, tf], u_init);
    elseif bc_method == "ipm"
        [t, u] = ode45(@sbp_ipm, [0, tf], u_init);
    else
        error('bc_method not recognized')
    end
end