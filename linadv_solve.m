% Simple Linear Advection Equation solver using SBP

% bc_method : boundary condition method ('sat', 'proj', 'ipm')
% n : grid size [0:1]
% tf : final time
% u_init : initial condition
% D1_func : SBP operator function
% u0_func : boundary condition at x=0, function of time
% u0_t_func : time derivative of boundary condition at x=0

function [t, u] = linadv_solve(odesolve, bc_method, n, tf, u_init, D1_func, u0_func, u0_t_func, sigma)
    [H, D1] = D1_func(n + 1);
    
    % scale SBP operator to grid size
    H = H / n;
    D1 = D1 * n;
    D1 = sparse(D1);
    
    % precalculate element 1, 1 of inv(H)
    invH = inv(H);
    invH_11 = invH(1, 1);

    %h = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
    h = 0.0001;
    if bc_method == "sat"
        [t, u] = odesolve(@(t, y) sbp_sat(t, y, D1, u0_func, invH_11), [0, tf], u_init, h);
    elseif bc_method == "proj"
        [t, u] = odesolve(@(t, y) sbp_proj(t, y, D1, u0_t_func), [0, tf], u_init, h);
    elseif bc_method == "ipm"
        [t, u] = odesolve(@(t, y) sbp_ipm(t, y, D1, u0_func, u0_t_func, sigma), [0, tf], u_init, h);
    else
        error('bc_method not recognized')
    end
end