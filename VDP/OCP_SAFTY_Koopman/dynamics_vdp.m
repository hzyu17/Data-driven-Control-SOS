function [F,G] = dynamics_vdp()

%Define the dynamics f and g for vdp
% pvar x1 x2
% x = [x1;x2];
% F = [x(2); (1 - x(1)^2) * x(2) - x(1)];
% G = [0; 1];

%Simple dynamics for double integrator multiple inputs
F = [0; 0];
G = [1, 0; 0, 1];

end

