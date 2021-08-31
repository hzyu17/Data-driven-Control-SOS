% Van der Pol Oscillator: d_x1 = x2; d_x2 = (1-x1^2)x2 - x1 + u
% USING SOSOPT!
% ---------------
% Main content:
% 1. SOS optimization problem formulation
% 2. SOS solving
% 3. result plotting, including vector field and several trajectories
% ---------------
% HISTORY
% 12/02 Tested results for the new polynomial generator function, 
%           still not perfect vector field;
% 12/02 Add comments and chapters to the code
% 12/02 Tested with L1, L2: divergent circle; with L1, L2-L1: strange
%           vector field. Need further testing.
% 12/03 Trying to find out what is still different from Hyungjin's code
%           besides gEDMD.
% 12/04 Test the SOSOPT and SOSTOOLS using the modified EDMD result.

% 01/12/2021 Test with gEDMD, with one single point and its derivatives on
% the trajectory, results is not satisfied. Realized no normalizing. Retried, same
% results.
% 01/13/2021 Try with two consecutive points on the trajectory. But what is
% the meaning of taking two points and their derivatives? Amswer: there is
% no need for two points. only one point + derivative is enough.

%01/13/2021 Find out what was wrong: the matrox apporxi of (div G) should
%be C_c' * (L2-L1) * Psi + (C_c' * Psi) .* [div(C_x' * (L2 - L1) * Psi)],
%not the same as EDMD, problem solved.


close all
clear all
clc

% flag choosing between EDMD and gEDMD
using_gEDMD = 1;

%% 
% -----------
% load data
% -----------
if using_gEDMD
%     load('gEDMD_res.mat')
    load('data/gEDMD_res_01_15_21_01_27_58.mat')
else
    load('sampling_data.mat')
    load('EDMD_res.mat')
end

%% 
% =======================
% Try with SOSOPT matlab toolbox 
% =======================

% ------------------------
% Seperated calculation
% ------------------------
% % ------------------------------------------------------------------
% % 12/03 Test with Hyungjin's L1 and L2 matrices: good result
% % -----------------------------------------------------------------
% load('external/VdP gEDMD/gEDMD/result_gEDMD.mat');
% Psi = Psi_poly;

%%
pvar x1 x2
x = [x1;x2];
Div_F = 0; Div_G = 0;
F = c_x'*L{1}*Psi;
if using_gEDMD
    G = c_x'*(L{2}-L{1})*Psi;
else %EDMD
    G = c_x'*L{2}*Psi;
end

for i1 = 1 : nx
    Div_Fi = diff(F(i1), x(i1));
    Div_F = Div_F + Div_Fi;
    Div_Gi = diff(G(i1), x(i1));
    Div_G = Div_G + Div_Gi;
end

Div_F.coefficient(find(abs(Div_F.coefficient) <= 1e-9)) = 0;
Div_G.coefficient(find(abs(Div_G.coefficient) <= 1e-9)) = 0;

% term1
term(1) = c_a'*L{1}*Psi + c_a'*Psi*Div_F;
% term2, term4
if using_gEDMD
    term(2) = c_c'*(L{2}-L{1})*Psi + c_c'*Psi*Div_G;
    term(4) = c_bc'*(L{2}-L{1})*Psi + c_bc'*Psi*Div_G;
else %EDMD
    term(2) = c_c'*L{2}*Psi + c_c'*Psi*Div_G;
    term(4) = c_bc'*L{2}*Psi + c_bc'*Psi*Div_G;
end
% term3
term(3) = c_ab'*L{1}*Psi + c_ab'*Psi*Div_F;

%%
% =================
% Formulate SOS constraint
% =================
% -------------------------------------
% Using polynomial representation
% -------------------------------------
constraint_sys = (1+alpha)*poly_b*sum(term(1:2)) - alpha*sum(term(3:4));
constraint_sys.coefficient(find(abs(constraint_sys.coefficient)<=1e-0)) = 0;

%%
% ----------------------------
% SOS constraint in solver
% ----------------------------
constraints = polyconstr; % initialize polynomial constraint
constraints(1) = constraint_sys >= 0; % Dual Lyapunov stability
constraints(2) = poly_a >= 0; % rho(x) >= 0

% -----------
% Solving
% -----------
opts = sosoptions;
[info, dopt, sossol] = sosopt(constraints, x, opts)

%%
% ============
% Get solutions 
% ============
% ax
C_a_Psi_sol = subs(c_a, dopt);
ax_sol = C_a_Psi_sol'*Psi;
% ax.coefficient = round(ax.coefficient,6);

% cx
C_c_Psi_sol = subs(c_c, dopt);
cx_sol = C_c_Psi_sol'*Psi
% cx.coefficient = round(cx.coefficient,6);

% abx
C_ab_Psi_sol = subs(c_ab, dopt);
ab_sol = C_ab_Psi_sol'*Psi;

% bcx
C_bc_Psi_sol = subs(c_bc, dopt);
bc_sol = C_bc_Psi_sol'*Psi;


%%
% ============
% Draw a vector field 
% ============
yl = -5.0;
yr = 5.0;
dy = 0.5;
[y1i, y2i] = meshgrid(yl:dy:yr, yl:dy:yr);
yi = [y1i(:), y2i(:)];
ucont = zeros(size(y1i));
vcont = zeros(size(y2i));

for i1 = 1:numel(y1i)
    xdot = feval(@(pt) control_dynamics_VDP(0, pt, x, cx_sol, ax_sol), yi(i1,: )');
    ucont(i1)= xdot(1);
    vcont(i1) = xdot(2);
end

figure
quiver(y1i, y2i, ucont, vcont);
xlabel('x_1');
ylabel('x_2');
axis tight equal;
if using_gEDMD
    title('vdp gEDMD')
else
    title('vdp EDMD')
end

% -------------------------
% plot result trajectories 
% -------------------------
pvar x1 x2
x = [x1; x2];
xinits = [-3,-3; -2,-3; 0,-3; 2,-3; 3,-3; -3,3; -2,3; 0,3; 2,3; 3,3];

% -----------------------------------------------------
% 2. solving ODE and stablizing using ux = cx / ax
% -----------------------------------------------------
tspan = [0 30];
for l1 = 1:length(xinits)
    xinit = xinits(l1, :)
    [tcont{l1}, ycont{l1}] = ode15s(@(t,pt) control_dynamics_VDP(t, pt, x, cx_sol, ax_sol), tspan, xinit');
end

hold on
grid on
hold on; cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont); hold off;

