function F = control_dynamics_VDP(t, pt, x, cx, ax, sin_x)
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

% VDP:
F = zeros(2,1);
F(1) = pt(2);
F(2) = (1 - pt(1)^2) * pt(2) - pt(1) + u;

% ======================
% Inverted pendulum:
% syms x1
% n = 3;
% sin_x1 = s2p(taylor(sin(x1), x1, 'Order', n));

% F = zeros(2,1);
% F(1) = pt(2);
% F(2) = 4*subs(sin_x, x, pt) + u;


% F = [pt(2); -pt(1) - pt(2) + (1/27)*pt(2)^3];
% G = [0; 0];

% ======================
% Boat:
% syms x1 x2
% n = 10;
% 
% sin_y = taylor(sin(0.5*x2), x2, 'Order', n);
% cos_x = taylor(cos(0.5*x1), x1, 'Order', n);
% 
% pvar x1 x2
% x = [x1;x2];
% 
% F_1 = subs(s2p(1 + 0.125*cos_x - 0.125* sin_y), x, pt);
% F = [F_1;  u];


end

