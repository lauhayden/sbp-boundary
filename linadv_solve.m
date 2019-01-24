% Simple Linear Advection Equation solver using SBP

% bc_method : boundary condition method ('sat', 'proj', 'ipm')
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

    if bc_method == "sat"
        [t, u] = ode45(@(t, y) sbp_sat(t, y, D1, u0_func, invH_11), [0, tf], u_init);
    elseif bc_method == "proj"
        [t, u] = ode45(@(t, y) sbp_proj(t, y, D1, u0_t_func), [0, tf], u_init);
    elseif bc_method == "ipm"
        [t, u] = ode45(@(t, y) sbp_ipm(t, y, D1, u0_func, u0_t_func, sigma), [0, tf], u_init);
    else
        error('bc_method not recognized')
    end
end