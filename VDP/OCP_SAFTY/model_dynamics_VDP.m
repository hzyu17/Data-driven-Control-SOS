function F = model_dynamics_VDP(t,x,u)
% dynamics of Van der Pol
if strcmp(class(x), 'polynomial')
    F = zeros(2, 1);
end
F(1, 1) = x(2);
F(2, 1) = (1 - x(1)^2) * x(2) - x(1) + u;
end

