function F = control_dynamics_VDP(t, pt, x, cx, ax)
% the stablizing control of Van Der Pol dynamics

rho = 28;
sigma = 10;
beta = 8/3;

u = subs(cx, x, pt) / subs(ax, x, pt);
F = zeros(3,1);
F(1) = sigma * (pt(2) - pt(1));
F(2) = pt(1) * (rho - pt(3)) - pt(2) + u;
F(3) = pt(1) * pt(2) - beta * pt(3);
end

