function F = model_dynamics_VDP(t,x,u)
% dynamics of Van der Pol

rho = 28;
sigma = 10;
beta = 8/3;

F = zeros(3, 1);
F(1) = sigma * (x(2) - x(1));
F(2) = x(1) * (rho - x(3)) - x(2) + u;
F(3) = x(1) * x(2) - beta * x(3);
end

