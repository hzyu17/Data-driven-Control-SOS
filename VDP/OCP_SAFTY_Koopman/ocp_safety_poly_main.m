function [outputArg1,outputArg2] = ocp_safety_rational_main(lambda, eps, geo)
%% SOS program creation or reading
load('sos_prog.mat')

%% calculate the OCP objective function in function of the c_w and c_a
% Psi_sym = p2s(Psi);

X_excld = geo.X_excld;
X_r_coord = geo.X_r_coord;
X_u_coord = geo.X_u_coord;

integral_vec_a_r = zeros(Qa, 1);

for li=1:Qa
    cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym);
    cost_fun2_hdl = matlabFunction(Psi_sym(li), 'Vars', [x1, x2]);
    integral_vec_a_r(li) = integral_vec_a_r(li) - integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_r_coord(4));
%     integral_vec_a_r(li) = integral_vec_a_r(li) - (integral2(cost_fun1_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun1_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4)));
    if lambda > 0
        if li == 1
            integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * (X_u_coord(2) - X_u_coord(1)) * (X_u_coord(4) - X_u_coord(3));
        else
            added = lambda * integral2(cost_fun1_hdl, X_u_coord(1), X_u_coord(2), X_u_coord(3), X_u_coord(4))
            integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * integral2(cost_fun2_hdl, X_u_coord(1), X_u_coord(2), X_u_coord(3), X_u_coord(4));
        end
    end
end
obj_sos = obj_sos + c_a_sym(1:Qa)' * integral_vec_a_r

integral_vec_w_r = zeros(Qw, 1);
integral_vec_w_r(1) = (X_0(2) - X_0(1)) * (X_0(4) - X_0(3));
for li=2:Qw
    cost_fun2_hdl = matlabFunction(Psi_sym(li), 'Vars', [x1, x2]);
    integral_vec_w_r(li) = integral_vec_w_r(li) - integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_r_coord(4));
%     integral_vec_w_r(li) = integral_vec_w_r(li) - (integral2(cost_fun2_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) +  integral2(cost_fun2_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4)));
end
obj_sos = obj_sos + c_w_sym(1:Qw)' * integral_vec_w_r;
sos_prog = sossetobj(sos_prog, obj_sos);


%%

syms x1 x2
x = [x1; x2];

term(1) = vec_div(F, x)*poly_a_sym + jacobian(poly_a_sym, x)*F;
term(2) = vec_div(G, x)*poly_c_sym + jacobian(poly_c_sym, x) * G;

constraint_sym = sum(term(1:2));
constraint_sym = vpa(constraint_sym, 9);
sos_prog = sosineq(sos_prog, constraint_sym);
sos_prog = sosineq(sos_prog, poly_a_sym);

c_s1 = mpvar('c_s1', [Qs1, 1]);
c_s1 = [c_s1(1:Qs1); zeros(total_deg-Qs1,1)];
s1 = p2s(c_s1' * Psi);
c_s1_sym = p2s(c_s1);
sos_prog = sosdecvar(sos_prog, c_s1_sym(1:Qs1));
sos_prog = sosineq(sos_prog, s1);

c_s3 = mpvar('c_s3', [Qs1, 1]);
c_s3 = [c_s3(1:Qs1); zeros(total_deg-Qs1,1)];
s3 = p2s(c_s3' * Psi);
c_s3_sym = p2s(c_s3);
sos_prog = sosdecvar(sos_prog, c_s3_sym(1:Qs1));
sos_prog = sosineq(sos_prog, s3);

% c_s4 = mpvar('c_s4', [Qs1, 1]);
% c_s4 = [c_s4(1:Qs1); zeros(total_deg-Qs1,1)];
% s4 = p2s(c_s4' * Psi);
% c_s4_sym = p2s(c_s4);
% sos_prog = sosdecvar(sos_prog, c_s4_sym(1:Qs1));
% sos_prog = sosineq(sos_prog, s4);
% 
% c_s5 = mpvar('c_s5', [Qs1, 1]);
% c_s5 = [c_s5(1:Qs1); zeros(total_deg-Qs1,1)];
% s5 = p2s(c_s5' * Psi);
% c_s5_sym = p2s(c_s5);
% sos_prog = sosdecvar(sos_prog, c_s5_sym(1:Qs1));
% sos_prog = sosineq(sos_prog, s5);

[X_0, X, X_d, X_init, X_u, X_r] = domain_definition();
% Constraints on X_init: div(pf+gu) > 1-eps
% ineq = constraint_sym - (1 - eps) - s1*X_init;
% ineq = vpa(ineq, 9);
% sos_prog = sosineq(sos_prog, ineq);

% constraints on X_u:div(pf+gu) < eps
% ineq3 = - constraint_sym + eps - s3*X_u;
% ineq3 = vpa(ineq3, 9);
% sos_prog = sosineq(sos_prog, ineq3);

% for polynomials, constraints on [X_d, X_d_out]: div(pf+gu)<eps)
% X_d_out = X_d + 2;
% ineq4 = -poly_a_sym + eps + s4*X_d - s5*X_d_out;
% ineq4 = vpa(ineq4, 9);
% sos_prog = sosineq(sos_prog, ineq4);

% c_s2 = mpvar('c_s2', [Qs1, 1]);
% c_s2 = [c_s2(1:Qs1); zeros(total_deg-Qs1,1)];
% s2 = p2s(c_s2' * Psi);
% c_s2_sym = p2s(c_s2);
% sos_prog = sosdecvar(sos_prog, c_s2_sym(1:Qs1));

% sos_prog = sosineq(sos_prog, s2);
% ineq2 = poly_a_sym - 10e4*(1 - eps) -  s2*X_r;
% ineq2 = vpa(ineq2, 9);
% sos_prog = sosineq(sos_prog, ineq2);

%% solve
disp('Solving sos program')

option.solver = 'sedumi';
option.params.tol = 1e-15;
[sos_prog, info] = sossolve(sos_prog, option);

disp('------- numerical issue -----------')
info.numerr

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

ax_sol_sym = vpa(ax_sol_sym, 9);
bx_sol_sym = vpa(bx_sol_sym, 9);
cx_sol_sym = vpa(cx_sol_sym, 9);
wx_sol_sym = vpa(wx_sol_sym, 9);

name = ['experiments/results_SOS_Poly_lbd_',num2str(lambda),'_eps_',num2str(eps),'.mat'];
save(name)

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


