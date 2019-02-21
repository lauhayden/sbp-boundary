function [t, u] = rk4_wrapper(odefun, tspan, y0, h)
    t = tspan(1):h:tspan(2);
    u = zeros(length(t), length(y0));
    u(1, :) = y0;
    for i=2:length(t)
        u(i, :) = rk4(h, odefun, t(i), u(i-1, :)');
    end
end