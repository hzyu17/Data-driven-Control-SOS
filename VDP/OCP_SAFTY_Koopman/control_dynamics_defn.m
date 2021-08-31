function F = control_dynamics_defn(dynamics_type, pt, x, cx, ax)
% the stablizing control of Van Der Pol dynamics: faster using polynomial
% subs
size_c = size(cx);
dim_m = size_c(1);
u = [];
for i_c = 1:dim_m
    if contains(class(cx(i_c)), 'polynomial')
%         subs(cx(i_c), x, pt) / subs(ax, x, pt);
        u = [u; subs(cx(i_c), x, pt) / subs(ax, x, pt)];
    elseif contains(class(cx), 'sym')
        syms x1 x2
        x = [x1; x2];
        f_hdl = matlabFunction(cx(i_c), 'Vars', x);
        feval(f_hdl, pt(1), pt(2))
        u = [u; feval(f_hdl, pt(1), pt(2))];
    end
end

dynamics = dynamics_definition('polynomial', dynamics_type);

F = double(subs(dynamics.F, x, pt));
G = double(dynamics.G);
F = F + G * double(u);

% % VDP:
% % F = zeros(2,1);
% % F(1) = pt(2);
% % F(2) = (1 - pt(1)^2) * pt(2) - pt(1) + u;
% 
% % Integrator
% F = [0;0];
% G = eye(dim_m);
% 
% % testing dynamics
% F = [-pt(1) - pt(2); pt(1) - pt(2)^3];
% G = [1, 0; 0, 1];
% % double(u)
% F = F + G * double(u);

% linear dynamics
% F = [-pt(1) ; - pt(2)];
% G = [1, 0; 0, 1];

% dynamics in Rantzer's paper
% F = [pt(2); -pt(1) + pt(1)^3./3 - pt(2)];
% G = [1, 0; 0, 1];

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

end

