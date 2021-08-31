% Van der Pol Oscillator: d_x1 = x2; d_x2 = (1-x1^2)x2 - x1 + u
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

clear all
close all
clc

%% 
% -----------
% load data
% -----------
load('sampling_data.mat')
load('EDMD_res.mat')

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
% -------------------------
% symbolic constraints
% -------------------------
% 1. SOS constraint
poly_c_sym = p2s(poly_c);

div_fa = c_a' * P{1} * Psi;
div_fa_sym = p2s(div_fa);

div_gc = c_c' * P{2} * Psi;
div_gc_sym = p2s(div_gc);

div_fab = c_ab' * P{1} * Psi;
div_fab_sym = p2s(div_fab);

div_gbc = c_bc' * P{2} * Psi;
div_gbc_sym = p2s(div_gbc);

b_sym = p2s(poly_b);

sos_poly_sym = (1+alpha) * b_sym * (div_fa_sym + div_gc_sym) - alpha * (div_fab_sym + div_gbc_sym);

%% 
% =============
% Use SOSTOOLS 
% =============

% --------------------------------------------
% Declare SOS program and independent 
% polynomial variables
% --------------------------------------------
syms x1 x2;
program2 = sosprogram([x1;x2]);

% ------------------------------------
% declare decision vairables: the c's
% ------------------------------------
c_c_sym = p2s(c_c);
program2 = sosdecvar(program2, [c_c_sym(1:Qc)]);

c_a_sym = p2s(c_a);
poly_a_sym = p2s(poly_a);
program2 = sosdecvar(program2, [c_a_sym(1:Qa)]);

% c_y = mpvar('y', [Qc, 1]);
% z = zeros(total_coeffs-Qc, 1);
% p = polynomial(z);
% c_y = [c_y;p];
% c_y_sym = p2s(c_y);
% poly_y_sym = p2s(c_y' * Phix);
% program2 = sosdecvar(program2, [c_y_sym(1:Qc)]);
% 
% c_s = mpvar('s', [Qa, 1]);
% z = zeros(total_coeffs-Qa, 1);
% p = polynomial(z);
% c_s = [c_s;p];
% c_s_sym = p2s(c_s);
% poly_s_sym = p2s(c_s' * Phix);
% program2 = sosdecvar(program2, [c_s_sym(1:Qa)]);

% -------------------
% SOS Constraints
% -------------------
constraint_poly = s2p(sos_poly_sym);
constraint_poly.coefficient(find(abs(constraint_poly.coefficient)<=1e-0)) = 0;
sos_poly_sym = p2s(constraint_poly);
program2 = sosineq(program2, sos_poly_sym);
program2 = sosineq(program2, poly_a_sym);

% ==========================
% constraints and the l1-norm objective
% ==========================
% % l1 norm for c(x)
% % for i=1:Qc
% % %     program2 = sosineq(program2, poly_y_sym - poly_c_sym);
% % %     program2 = sosineq(program2, poly_y_sym + poly_c_sym);
% %     program2 = sosineq(program2, c_y_sym(i) + c_c_sym(i));
% %     program2 = sosineq(program2, c_y_sym(i) - c_c_sym(i));
% % end
% % 
% % % =========
% % l1 norm for a(x)
% % % =========
% % for i=1:Qa
% % %     program2 = sosineq(program2, poly_s_sym - poly_c_sym);
% % %     program2 = sosineq(program2, poly_s_sym + poly_c_sym);
% %     program2 = sosineq(program2, c_s_sym(i) + c_a_sym(i));
% %     program2 = sosineq(program2, c_s_sym(i) - c_a_sym(i));
% % end
% 
% % ====================
% % l1-norm Objective function 
% % ====================
% % c_obj = ones(Qc+Qa, 1);
% % objective1 = c_obj(1:Qc)' * c_y_sym(1:Qc) + c_obj(Qc+1:Qc+Qa)' * c_s_sym(1:Qa)
% % program2 = sossetobj(program2, objective1);
% 
% ============
% Solve the program 
% ============
[program2, info] = sossolve(program2);

% =========== 
% Get SOS solution 
% ===========
%  -----------------------------------
% Solution of decision variables: c_c, c_y, c_a
% ------------------------------------
% c_c_sym_sol = sosgetsol(program2, c_c_sym)
% c_y_sym_sol = sosgetsol(program2, c_y_sym)
% c_a_sym_sol = sosgetsol(program2, c_a_sym)
sol_poly_main_sym = sosgetsol(program2, sos_poly_sym);
% 
% % -------------------------
% % L1 objective solution
% % -------------------------
% % objective1 = sosgetsol(program2, objective1)
% 
% -------------------------------
% polynomial solution: cx, ax.
% -------------------------------
cx_sol_sym = sosgetsol(program2, poly_c_sym)
cx_sol = s2p(cx_sol_sym);
ax_sol_sym = sosgetsol(program2, poly_a_sym);
ax_sol = s2p(ax_sol_sym);

pvar x1 x2;
x = [x1;x2];

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
tspan = [0 20];
for l1 = 1:length(xinits)
    xinit = xinits(l1, :)
    [tcont{l1}, ycont{l1}] = ode15s(@(t,pt) control_dynamics_VDP(t, pt, x, cx_sol, ax_sol), tspan, xinit');
end

hold on
grid on
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont); hold off;
