function [F,G] = dynamics_vdp()

%Define the dynamics f and g for vdp
pvar x1 x2
x = [x1;x2];
F = [x(2); (1 - x(1)^2) * x(2) - x(1)];
G = [0; 1];

% dynamics for planar quadrotor
% m = 1;
% I = 1;
% g = 9.81;
% r = 2;
% syms x1 x2 x3 x4 x5 x6 
% x = [x1;x2;x3;x4;x5;x6];
% F = [x4; 
%         x5;
%         x6;
%         0;
%         -g;
%         0];
% G = [0,0; 0,0; 0,0;
%         -sin(x3)/m, -sin(x3)/m;
%         cos(x3)/m, cos(x3)/m;
%         r/I, -r/I];

% F = [x(2)-x(1)^3+x(1)^2; 0];
% G = [0; 1];
end

