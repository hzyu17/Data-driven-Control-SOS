function [F,G] = dynamics_craft()

%Define the dynamics f and g for vdp
pvar x1 x2 x3
x = [x1;x2;x3];

% Dynamics for a rigid rotating spacecraft
J_1 = 2;
J_2 = 3;

F = [(J_2 -  J_3) / J_1 * x2 * x3;
    (J_3 - J_1)/J_2 * x3 * x1;
    (J_1 - J_2)/J_3 * x1 * x2];

G = [1/J_1, 0, 0;
    0, 1/J_2, 0;
    0, 0, 1/J_3];

end

