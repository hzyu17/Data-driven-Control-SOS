function [out] = create_sosprog(poly)
%%
% polynomials a, b, c, q, R, w: if it does not exist, please first run polynomials_definition.m
% close all 
% clear all
% clc

load('polynomials_def.mat')

% domains
%%
geo = domain_definition();
X = geo.X;
X = geo.X;
X_d = geo.X_d;
X_init = geo.X_init;
X_r = geo.X_r;
X_u = geo.X_u;

%%
% dynamics
[F, G] = dynamics_vdp();

F = p2s(F);
F = vpa(F, 9);

syms x1 x2
x = [x1; x2];

%% sos constraints on OCP equivalent matrix SDP
% add poly_w as decision variables
syms x1 x2
sos_prog = sosprogram([x1, x2]);

c_c_sym = p2s(c_c);
sos_prog = sosdecvar(sos_prog, [c_c_sym(1:Qc)]);

% c_c1_sym = p2s(c_c1);
% sos_prog = sosdecvar(sos_prog, [c_c1_sym(1:Qc1)]);

c_a_sym = p2s(c_a);
% poly_a_sym = p2s(poly_a);
sos_prog = sosdecvar(sos_prog, [c_a_sym(1:Qa)]);

% sos constraints on density and a(x)
% sos_prog = sosineq(sos_prog, constraint_sym);
% sos_prog = sosineq(sos_prog, poly_a_sym);

c_w_sym = p2s(c_w);
sos_prog = sosdecvar(sos_prog, [c_w_sym(1:Qw)]);

% PSD Matrix M of polynomials
size_M = 1+length(poly_c);
M = sym('M', [2,2]);
M(1,1) = poly_w_sym;
M(1,2) = poly_c_sym;
M(2,1) = poly_c_sym;
M(2:end, 2:end) = R^-1 * poly_a_sym; % Should be inv(R) in general cases.

sos_prog = sosmatrixineq(sos_prog, M, 'quadraticMineq');

%% Base integral
Psi_sym = p2s(Psi);
integral_vec_a = zeros(Qa, 1);

X_init = geo.X_init;
X_excld = geo.X_excld;
X_r_coord = geo.X_r_coord;
X_u_coord = geo.X_u_coord;
X_0 = geo.X_0;

for li=1:Qa
    if poly == 1
        cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym);
        integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_0(4));
    else % rational
        cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym / (poly_b_sym^Alph));
        integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun1_hdl, X_excld(2), X_0(2), X_0(3), X_0(4)) + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
%         integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun1_hdl, X_excld(2), X_0(2), X_0(3), X_0(4));
%                              + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    end
end
obj_sos = c_a_sym(1:Qa)' * integral_vec_a;

integral_vec_w = zeros(Qw, 1);
integral_vec_w(1) = (X_0(2) - X_0(1)) * (X_0(4) - X_0(3));
for li=2:Qw
    if poly == 1
        cost_fun2_hdl = matlabFunction(Psi_sym(li), 'Vars', [x1, x2]);
        integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_0(4));
    else % rational
        cost_fun2_hdl = matlabFunction(Psi_sym(li) / (poly_b_sym^Alph));
        integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun2_hdl, X_excld(2), X_0(2), X_0(3), X_0(4)) + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
%         integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun2_hdl, X_excld(2), X_0(2), X_0(3), X_0(4));
%                              + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    end
end
obj_sos = obj_sos + c_w_sym(1:Qw)' * integral_vec_w;


%% sos program for safety verification
% define a new sos programsyms x1 x2
% syms x1 x2
% sos_safety_prog = sosprogram([x1, x2]);
% 
% deg_safe_rho = 12;
% Q_safe_rho = nchoosek(deg_safe_rho+nx, nx);
% c_safe_rho = mpvar('c_safe_rho', [Q_safe_rho, 1]);
% c_safe_rho = [c_safe_rho(1:Q_safe_rho); zeros(nPsi-Q_safe_rho, 1)];
% c_safe_rho_sym = p2s(c_safe_rho);
% poly_safe_rho_sym = p2s(c_safe_rho' * Psi);
% 
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_rho_sym(1:Q_safe_rho));
% 
% deg_safe_p1 = polynomialDegree(X);
% deg_safe_s1 = deg_safe_rho - deg_safe_p1;
% Q_safe_s1 = nchoosek(deg_safe_s1+nx, nx);
% c_safe_s1 = mpvar('c_safe_s1', [Q_safe_s1, 1]);
% c_safe_s1 = [c_safe_s1(1:Q_safe_s1); zeros(nPsi-Q_safe_s1,1)];
% poly_safe_s1 = p2s(c_safe_s1' * Psi);
% c_safe_s1_sym = p2s(c_safe_s1);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s1_sym(1:Q_safe_s1));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s1);
% 
% deg_safe_p2 = polynomialDegree(X);
% deg_diff =  total_deg - (deg_safe_rho + 1 + 2*deg_a - deg_safe_p2)
% deg_safe_s2 = total_deg - deg_safe_p2;
% Q_safe_s2 = nchoosek(deg_safe_s2+nx, nx);
% c_safe_s2 = mpvar('c_safe_s2', [Q_safe_s2, 1]);
% c_safe_s2 = [c_safe_s2(1:Q_safe_s2); zeros(nPsi-Q_safe_s2, 1)];
% poly_safe_s2 = p2s(c_safe_s2' * Psi);
% c_safe_s2_sym = p2s(c_safe_s2);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s2_sym(1:Q_safe_s2));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s2);
% 
% deg_safe_p3 = polynomialDegree(X_r);
% deg_safe_s3 = total_deg - deg_safe_p3;
% Q_safe_s3 = nchoosek(deg_safe_s3+nx, nx);
% c_safe_s3 = mpvar('c_safe_s3', [Q_safe_s3, 1]);
% c_safe_s3 = [c_safe_s3(1:Q_safe_s3); zeros(nPsi-Q_safe_s3, 1)];
% poly_safe_s3 = p2s(c_safe_s3' * Psi);
% c_safe_s3_sym = p2s(c_safe_s3);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s3_sym(1:Q_safe_s3));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s3);
% 
% deg_safe_p4 = polynomialDegree(X_d);
% deg_safe_s4 = deg_safe_rho - deg_safe_p4;
% Q_safe_s4 = nchoosek(deg_safe_s4+nx, nx);
% c_safe_s4 = mpvar('c_safe_s4', [Q_safe_s4, 1]);
% c_safe_s4 = [c_safe_s4(1:Q_safe_s4); zeros(nPsi-Q_safe_s4, 1)];
% poly_safe_s4 = p2s(c_safe_s4' * Psi);
% c_safe_s4_sym = p2s(c_safe_s4);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s4_sym(1:Q_safe_s4));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s4);
% 
% deg_safe_p5 = polynomialDegree(X);
% deg_safe_s5 = deg_safe_rho - deg_safe_p5;
% Q_safe_s5 = nchoosek(deg_safe_s5+nx, nx);
% c_safe_s5 = mpvar('c_safe_s5', [Q_safe_s5, 1]);
% c_safe_s5 = [c_safe_s5(1:Q_safe_s5); zeros(nPsi-Q_safe_s5, 1)];
% poly_safe_s5 = p2s(c_safe_s5' * Psi);
% c_safe_s5_sym = p2s(c_safe_s5);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s5_sym(1:Q_safe_s5));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s5);
% 
% deg_safe_p6 = polynomialDegree(X_u);
% deg_safe_s6 = deg_safe_rho - deg_safe_p6;
% Q_safe_s6 = nchoosek(deg_safe_s6+nx, nx);
% c_safe_s6 = mpvar('c_safe_s6', [Q_safe_s6, 1]);
% c_safe_s6 = [c_safe_s6(1:Q_safe_s6); zeros(nPsi-Q_safe_s6, 1)];
% poly_safe_s6 = p2s(c_safe_s6' * Psi);
% c_safe_s6_sym = p2s(c_safe_s6);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s6_sym(1:Q_safe_s6));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s6);
% 
% deg_safe_p7 = polynomialDegree(X_r);
% deg_safe_s7 = deg_safe_rho - deg_safe_p7;
% Q_safe_s7 = nchoosek(deg_safe_s7+nx, nx);
% c_safe_s7 = mpvar('c_safe_s7', [Q_safe_s7, 1]);
% c_safe_s7 = [c_safe_s7(1:Q_safe_s7); zeros(nPsi-Q_safe_s7, 1)];
% poly_safe_s7 = p2s(c_safe_s7' * Psi);
% c_safe_s7_sym = p2s(c_safe_s7);
% sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s7_sym(1:Q_safe_s7));
% sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s7);
load('safety_prog.mat')
save('sos_prog')
end

