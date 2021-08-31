function F = control_dynamics_VDP(t, pt, x, cx, ax)
% the stablizing control of Van Der Pol dynamics: faster using polynomial
% subs
if contains(class(cx), 'polynomial')
    u = subs(cx, x, pt) / subs(ax, x, pt);
elseif contains(class(cx), 'sym')
    syms x1 x2
    x = [x1; x2];
    u_sym = vpa(cx/ax, 9);
    f_hdl = matlabFunction(u_sym, 'Vars', x);
    u = feval(f_hdl, pt(1), pt(2));
end
F = zeros(2,1);
F(1) = pt(2);
F(2) = (1 - pt(1)^2) * pt(2) - pt(1) + u;
end

