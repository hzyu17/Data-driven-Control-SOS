%% Optimal control for nonlinear dynamical systems with L2 cost and safty constraints using SOS
% assume we know the dynamics
% process: 
% OCP Part:
% 1. Read the dynamics f(x) from functions or Koopman
% 2. Read polynomials and Create SOS program
% 2. Define the reward matrix q(x), R(x), and calculate integral OCP cost
% 3. Define the Matrix equality and SOS feasibility problem
% 4. Add the density SOS constraints

% Safty Part:
% 1. Read the initial sets, unsafe sets, and reachable sets 
% 2. Formulate the SOS constraints w.r.t. X_0, X_u and X_r
% 3. Solving the whole program

% Log
% 02/20/2021: Start from OCP problem: recover the polynomials without
% Koopman
% 02/21/2021: Formulate and solve the OCP problem in matrix SDP in SOSTOOLS.
% 03/07/2021: Formulate the safety cost in the formulation

% close all 
clear all
clc
load('polynomials_definition.mat')

% domains
%%
[X_0, X, X_d, X_init, X_u, X_r] = domain_definition();
% X = 25 - x1^2 - x2^2;
% X_d = 36 - x1^2 - x2^2;
% X_init = 4 - (x1+2)^2 - (x2-2)^2;
% X_r = 4 - (x1-2)^2 - (x2+2)^2;
% X_u = 1 - (x1-1)^2 - (x2-1)^2;

%%
% dynamics
[F, G] = dynamics_vdp();
poly_a_sym = p2s(poly_a);
poly_b_sym = p2s(poly_b);
poly_c_sym = p2s(poly_c);

poly_w_sym = p2s(poly_w);

F = p2s(F);
G = p2s(G);
F = vpa(F, 9);
G = vpa(G, 9);
syms x1 x2
x = [x1; x2];
% constraint_sym = densityConstraint(poly_a_sym, poly_b_sym, poly_c_sym, Alph, F, G, nx);

constraint_sym = poly_b_sym*sum(diag(jacobian(F*poly_a_sym+G*poly_c_sym,x))) - Alph*jacobian(poly_b_sym,x)*(poly_a_sym*F+G*poly_c_sym);

% Rho = poly_a_sym / ((poly_b_sym)^Alph);
% Rho_bar = poly_c_sym / ((poly_b_sym)^Alph);
% LHS = Rho*F + Rho_bar*G;
% constraint_sym1 = diff(LHS(1), x1) + diff(LHS(2), x2);

% With model:
% rho=a/b^Alph, rho_bar = c/b^Alph

% Rho = poly_a_sym / ((poly_b_sym)^Alph);
% Rho_bar = poly_c_sym / ((poly_b_sym)^Alph);
% LHS = Rho*F + Rho_bar*G;
% constraint_sym = diff(LHS, x1) + diff(LHS, x2);

%% sos constraints on OCP equivalent matrix SDP
% add poly_w as decision variables
syms x1 x2
sos_prog = sosprogram([x1, x2]);
c_c_sym = p2s(c_c);
sos_prog = sosdecvar(sos_prog, [c_c_sym(1:Qc)]);

c_a_sym = p2s(c_a);
poly_a_sym = p2s(poly_a);
sos_prog = sosdecvar(sos_prog, [c_a_sym(1:Qa)]);

% sos constraints on density and a(x)
% sos_prog = sosineq(sos_prog, constraint_sym);
sos_prog = sosineq(sos_prog, poly_a_sym);

c_w_sym = p2s(c_w);
sos_prog = sosdecvar(sos_prog, [c_w_sym(1:Qw)]);

% PSD Matrix M of polynomials
size_M = 1+length(poly_c);
M = sym('M', [2,2]);
M(1,1) = poly_w_sym;
M(1,2:end) = poly_c_sym';
M(2:end,1) = poly_c_sym;
M(2:end, 2:end) = R^-1 * poly_a_sym; % Should be inv(R) in general cases.

sos_prog = sosmatrixineq(sos_prog, M, 'quadraticMineq');

%% SOS program creation
% create_sosprog();
% load('sos_prog.mat')

%%
% calculate the OCP objective function in function of the c_w and c_a
Psi_sym = p2s(Psi);
integral_vec_a = zeros(Qa, 1);
X_excld = [-0.1, 0.1, -0.1, 0.1];
X_r_coord = [0, 2, -2, 0];
X_u_coord = [0 1 0 1];
lambda = 1;
for li=1:Qa
    cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym / (poly_b_sym^Alph));
    integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun1_hdl, X_excld(2), X_0(2), X_0(3), X_0(4))
                             + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    integral_vec_a(li) = integral_vec_a(li) - integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_r_coord(4))
    integral_vec_a(li) = integral_vec_a(li) + integral2(cost_fun1_hdl, X_u_coord(1), X_u_coord(2), X_u_coord(3), X_u_coord(4))
end
obj_sos = c_a_sym(1:Qa)' * integral_vec_a;

integral_vec_w = zeros(Qw, 1);
for li=1:Qw
    cost_fun2_hdl = matlabFunction(Psi_sym(li) / (poly_b_sym^Alph));
    integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun2_hdl, X_excld(2), X_0(2), X_0(3), X_0(4))
                             + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    integral_vec_w(li) = integral_vec_w(li) - integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_r_coord(4))
end
obj_sos = obj_sos + c_w_sym(1:Qw)' * integral_vec_w;
sos_prog = sossetobj(sos_prog, obj_sos);

%% Add equality constraint: div(f*rho + g*rho_bar) = h(x)
% LHS = Rho*F + Rho_bar*G;
% div_LHS = diff(LHS, x1) + diff(LHS, x2);
eps = 1e-9;

total_deg = nchoosek(highest_order+nx, nx);
deg_p1 = polynomialDegree(X);
deg_s1 = highest_order - deg_p1;
Qs1 = nchoosek(deg_s1+nx, nx);

%%
% sos_prog = sosineq(sos_prog, constraint_sym);

% syms x1 x2
% x = [x1; x2];
% constraint_sym1 = poly_b_sym*sum(diag(jacobian(F*poly_a_sym+G*poly_c_sym,x))) - Alph*jacobian(poly_b_sym,x)*(poly_a_sym*F+G*poly_c_sym);

% sos_prog = sosineq(sos_prog, constraint_sym);

% Constraints on X_init: div(pf+gu) > 1-eps
c_s1 = mpvar('c_s1', [Qs1, 1]);
c_s1 = [c_s1(1:Qs1); zeros(total_deg-Qs1,1)];
s1 = p2s(c_s1' * Psi);
c_s1_sym = p2s(c_s1);
sos_prog = sosdecvar(sos_prog, c_s1_sym(1:Qs1));
sos_prog = sosineq(sos_prog, s1);

ineq = constraint_sym - (1 - eps) * poly_b_sym^(Alph+1) - s1*X_init;
ineq = vpa(ineq, 9);
sos_prog = sosineq(sos_prog, ineq);

% c_s2 = mpvar('c_s2', [Qs1, 1]);
% c_s2 = [c_s2(1:Qs1); zeros(total_deg-Qs1,1)];
% s2 = p2s(c_s2' * Psi);
% c_s2_sym = p2s(c_s2);
% sos_prog = sosdecvar(sos_prog, c_s2_sym(1:Qs1));

% sos_prog = sosineq(sos_prog, s2);

% ineq2 = -constraint_sym + (1 + eps) * poly_b_sym^(Alph+1) - s2*X_init;
% ineq2 = vpa(ineq2, 9);
% sos_prog = sosineq(sos_prog, ineq2);

% constraints on X_u:div(pf+gu) < eps
c_s3 = mpvar('c_s3', [Qs1, 1]);
c_s3 = [c_s3(1:Qs1); zeros(total_deg-Qs1,1)];
s3 = p2s(c_s3' * Psi);
c_s3_sym = p2s(c_s3);
sos_prog = sosdecvar(sos_prog, c_s3_sym(1:Qs1));
sos_prog = sosineq(sos_prog, s3);

ineq3 = -constraint_sym + eps * poly_b_sym^(Alph+1) + s3*X_init;
ineq3 = vpa(ineq3, 9);
sos_prog = sosineq(sos_prog, ineq3);

%%
% solve
[sos_prog, info] = sossolve(sos_prog);

%%
% get solution: cx, ax.
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym)
cx_sol = s2p(cx_sol_sym);
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym)
ax_sol = s2p(ax_sol_sym);
wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
wx_sol = s2p(wx_sol_sym);
bx_sol_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol = s2p(bx_sol_sym);

% ax_sol_sym = vpa(ax_sol_sym, 5);
% bx_sol_sym = vpa(bx_sol_sym, 5);
% cx_sol_sym = vpa(cx_sol_sym, 5);
% wx_sol_sym = vpa(wx_sol_sym, 5);

save('results_SOS')

%% add the safety verification to SOS program
% rho=a/b^Alph, rho_bar = c/b^Alph

% epsilon = 1e-6;
% total_deg = nchoosek(highest_order+nx, nx);
% 
% deg_p1 = polynomialDegree(X);
% deg_s1 = highest_order - deg_p1;
% Qs1 = nchoosek(deg_s1+nx, nx);
% c_s1 = mpvar('c_s1', [Qs1, 1]);
% c_s1 = [c_s1(1:Qs1); zeros(total_deg-Qs1,1)];
% c_s1_sym = p2s(c_s1);
% sos_prog = sosdecvar(sos_prog, c_s1_sym(1:Qs1));
% 
% deg_p2 = polynomialDegree(X_u);
% deg_s2 = highest_order - deg_p2;
% Qs2 = nchoosek(deg_s2+nx, nx);
% c_s2 = mpvar('c_s2', [Qs2, 1]);
% c_s2 = [c_s2(1:Qs2); zeros(total_deg-Qs2, 1)];
% c_s2_sym = p2s(c_s2);
% sos_prog = sosdecvar(sos_prog, c_s2_sym(1:Qs2));
% 
% deg_p3 = polynomialDegree(X_r);
% deg_s3 = highest_order - deg_p3;
% Qs3 = nchoosek(deg_s3+nx, nx);
% c_s3 = mpvar('c_s3', [Qs3, 1]);
% c_s3 = [c_s3(1:Qs3); zeros(total_deg-Qs3, 1)];
% c_s3_sym = p2s(c_s3);
% sos_prog = sosdecvar(sos_prog, c_s3_sym(1:Qs3));

%% Drawings

plot_result()
% plot_result_safety(X_0, X_init, X_u, X_r)



