function F = control_dynamics_VDP(t, pt, x, cx, ax)
% the stablizing control of Van Der Pol dynamics
u = subs(cx, x, pt) / subs(ax, x, pt);
F = zeros(2,1);
F(1) = pt(2);
F(2) = (1 - pt(1)^2) * pt(2) - pt(1) + u;
end

