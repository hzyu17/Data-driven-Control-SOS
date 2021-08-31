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

% 01/17/2021 Do some test on the scalar inequality of the SOSTOOL to
% realize the coefficient l1-norm constraints. Findings: the l1-norm
% variable c_c_abs is very close to c_c, and close to zero. Although in
% some dimensions c_c_abs - c_c < 0. Is this result reasonable? Do a test!

% 01/17/2021 Write another test script to test the l1-norm objective optimization.

clc
clear all
close all

% flag choosing between EDMD and gEDMD
using_gEDMD = 1;


%% 
% -----------
% load data
% -----------

if using_gEDMD
%     load('gEDMD_res.mat')
    load('data/gEDMD_res_01_15_21_01_27_58.mat') %bx = LQR
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
% =================
% Formulate SOS constraint
% =================
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

% -------------------------------------
% polynomial representation
% -------------------------------------
poly_a_sym = p2s(poly_a);
poly_b_sym = p2s(poly_b);
poly_c_sym = p2s(poly_c);
constraint_poly = (1+alpha)*poly_b*sum(term(1:2)) - alpha*sum(term(3:4));
constraint_poly.coefficient(find(abs(constraint_poly.coefficient)<=1e-0)) = 0;
constraint_sym = p2s(constraint_poly);

%% 
% =============
% SOSTOOLS 
% =============

% --------------------------------------------------------
% Create SOS program and independent variables
% --------------------------------------------------------
syms x1 x2;
program2 = sosprogram([x1;x2]);

% -----------------------------------------------
% declare decision vairables: the c's and a's
% -----------------------------------------------

% ======================
% density SOS vars and constraints
% ======================
c_c_sym = p2s(c_c);
program2 = sosdecvar(program2, [c_c_sym(1:Qc)]);

c_a_sym = p2s(c_a);
poly_a_sym = p2s(poly_a);
program2 = sosdecvar(program2, [c_a_sym(1:Qa)]);

program2 = sosineq(program2, constraint_sym);
program2 = sosineq(program2, poly_a_sym);

% =================================
% l1-norm related variables, objective and constraints
% =================================
total_coeffs = size(Psi, 1);

c_c_abs = mpvar('y', [Qc, 1]);
z = zeros(total_coeffs-Qc, 1);
p = polynomial(z);
c_c_abs = [c_c_abs;p];
c_c_abs_sym = p2s(c_c_abs);
poly_c_abs_sym = p2s(c_c_abs' * Psi);
program2 = sosdecvar(program2, [c_c_abs_sym(1:Qc)]);

c_a_abs = mpvar('s', [Qa, 1]);
z = zeros(total_coeffs-Qa, 1);
p = polynomial(z);
c_a_abs = [c_a_abs;p];
c_a_abs_sym = p2s(c_a_abs);
poly_a_abs_sym = p2s(c_a_abs' * Psi);
program2 = sosdecvar(program2, [c_a_abs_sym(1:Qa)]);

% l1 norm for C_c and C_a
for i=1:Qc
    program2 = sosineq(program2, c_c_abs_sym(i) + c_c_sym(i));
    program2 = sosineq(program2, c_c_abs_sym(i) - c_c_sym(i));
end
for i=1:Qa
    program2 = sosineq(program2, c_a_abs_sym(i) + c_a_sym(i));
    program2 = sosineq(program2, c_a_abs_sym(i) - c_a_sym(i));
end

% L1 - sparsity objectives: ||a||_1 + ||c||_1
c_obj = ones(Qc+Qa, 1);
objective1 = c_obj(1:Qc)' * c_c_abs_sym(1:Qc) + c_obj(Qc+1:Qc+Qa)' * c_a_abs_sym(1:Qa)
program2 = sossetobj(program2, objective1);

% ============
% Solve the program 
% ============
[program2, info] = sossolve(program2);

% =========== 
% Get SOS solution 
% ===========
%  ------------------------------------------
% Solution of decision variables: c_c, c_a
% -------------------------------------------
c_c_sym_sol = sosgetsol(program2, c_c_sym)
c_a_sym_sol = sosgetsol(program2, c_a_sym)
sol_poly_main_sym = sosgetsol(program2, constraint_sym);

%  ---------------------------------------------------
% Solution of l1-norm decision variables: c_c_abs
% ----------------------------------------------------
c_c_abs_sym_sol = sosgetsol(program2, c_c_abs_sym)

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
