function F = model_dynamics_VDP(t,x,u)
% dynamics of Van der Pol
if strcmp(class(x), 'polynomial')
    F = zeros(2, 1);
else
F(1, 1) = x(2);
F(2, 1) = sin(x(1)) - 0.5*x(2) + u;
end

