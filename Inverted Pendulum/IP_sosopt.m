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
    load('data/gEDMD_res_01_15_21_01_29_46.mat')
else
    load('sampling_data.mat')
    load('EDMD_res.mat')
end

%%
% ---------------
% Test matrix K
% ---------------
% c_x' * K * Phix = [x1+v*dt; x2+v*dt] 
% for Van der Pol dynamcs, c_x'*K*Psi should be:
% c_x' * K{1} * Phix
% [x1+x2*dt; x2+[(1-x1^2)x2 - x1]*dt] = [x1+0.01x2; 1.01x2-0.01x1-0.01x1^2*x2]
% c_x' * K{2} * Phix
% [x1+0.01x2; 1.01x2-0.01x1-0.01x1^2*x2+0.01]


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
else
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
% term2
if using_gEDMD
    term(2) = c_c'*(L{2}-L{1})*Psi + c_c'*Psi*Div_G;
else %EDMD
    term(2) = c_c'*L{2}*Psi + c_c'*Psi*Div_G;
end
% term3
term(3) = c_ab'*L{1}*Psi + c_ab'*Psi*Div_F;
% term4
if using_gEDMD
    term(4) = c_bc'*(L{2}-L{1})*Psi + c_bc'*Psi*Div_G;
else %EDMD
    term(4) = c_bc'*L{2}*Psi + c_bc'*Psi*Div_G;
end

%%
% -------------------------
% One-step calculation
% -------------------------
% P{1} = L{1} + (div1.* eye(length(Psi)));
% P{2} = L{2} + (div2.* eye(length(Psi)));
%
% term(1) = c_a'*P{1}*Phix;
% % term1_sym = p2s(term(1));
% term(2) = c_c'*P{2}*Phix;
% % term2_sym = p2s(term(2));
% term(3) = c_ab'*P{1}*Phix;
% % term3_sym = p2s(term(3));
% term(4) = c_bc'*P{2}*Phix;
% % term4_sym = p2s(term(4));

%%
% --------------------------------------------
% Test calculation, numerical imprecision
% --------------------------------------------
% c_ab' * (Div_F.* eye(total_coeffs)) * Phix - c_ab' * Phix * Div_F;
% aa = c_ab' * P{1} * Phix - (c_ab' * L{1} * Phix + c_ab' * (Div_F.* eye(total_coeffs)) * Phix)
% bb = c_ab' * P{1} * Phix - c_ab' * L{1} * Phix - c_ab' * (Div_F.* eye(total_coeffs)) * Phix
% 
% c_bc'*P{2}*Phix - (c_bc' * L{2} * Phix + c_bc' * Phix * div2);


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
    xdot = feval(@(pt) control_dynamics_IP(0, pt, x, cx_sol, ax_sol), yi(i1,: )');
    ucont(i1)= xdot(1);
    vcont(i1) = xdot(2);
end

figure
quiver(y1i, y2i, ucont, vcont);
xlabel('x_1');
ylabel('x_2');
axis tight equal;
if using_gEDMD
    title('Inverted Pendulum gEDMD')
else
    title('Inverted Pendulum EDMD')
end

% -------------------------
% plot result trajectories 
% -------------------------
pvar x1 x2
x = [x1; x2];
xinits = [-3,-3; -2,-3; 0,-3; 2,-3; 3,-3; -3,3; -2,3; 0,3; 2,3; 3,3];

% -----------------------------------------------------
% 1. Plotting linearized dynamics using ux = cx / ax
% -----------------------------------------------------
% figure
% hold on
% grid on
% x_pts = [];
% xdot = [x2; (1-x1^2)*x2 - x1];
% for kk = 1:4
%     for i=1:4000
%         if i == 1
%             x_pts = xinits(:, kk);
%         else
%             u_solved_value = subs(cx_sol, x, x_pts(1:2, i-1)) / subs(ax_sol, x, x_pts(1:2, i-1))
%             x_dot_poly_value = subs(xdot, x, x_pts(1:2, i-1));
%             x_dot_poly_value(2) = x_dot_poly_value(2) + u_solved_value;
%             x_dot_pt = x_dot_poly_value;
%             x_pts(1:2, i) = x_pts(1:2, i-1) + x_dot_pt.*dt;
%             plot(x_pts(1,i), x_pts(2,i),'.')
%         end
%     end
% end

% -----------------------------------------------------
% 2. solving ODE and stablizing using ux = cx / ax
% -----------------------------------------------------
tspan = [0 30];
for l1 = 1:length(xinits)
    xinit = xinits(l1, :)
    [tcont{l1}, ycont{l1}] = ode15s(@(t,pt) control_dynamics_IP(t, pt, x, cx_sol, ax_sol), tspan, xinit');
end

hold on
grid on
hold on; cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont); hold off;
% xlim([ys ye]); ylim([ys ye]);
% ylabel('$x_2$','Interpreter','Latex','FontSize', 38);
% xlabel('$x_1$','Interpreter','Latex','FontSize', 38);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 38, 'TickLabelInterpreter', 'Latex');
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 38, 'TickLabelInterpreter', 'Latex');
% savefig('result_VdP.fig');
