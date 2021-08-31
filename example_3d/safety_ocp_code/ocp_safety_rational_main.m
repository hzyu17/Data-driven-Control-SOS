function [outputArg1,outputArg2] = ocp_safety_rational_main(lambda, eps, geo)
%% SOS program creation or reading
disp('Loading predefined sos program')
load('sos_prog.mat')

%% calculate the OCP objective function in function of the c_w and c_a
% Psi_sym = p2s(Psi);

X_excld = geo.X_excld;
X_r_coord = geo.X_r_coord;
X_u_coord = geo.X_u_coord;

integral_vec_a_r = zeros(Qa, 1);

for li=1:Qa
    cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym / (poly_b_sym^Alph));
    cost_fun2_hdl = matlabFunction(Psi_sym(li) / (poly_b_sym^Alph), 'Vars', [x1, x2]);
    integral_vec_a_r(li) = integral_vec_a_r(li) - integral2(cost_fun1_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) + ... 
                                                                  integral2(cost_fun1_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + ...
                                                                  integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + ...
                                                                  integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4));
%     integral_vec_a_r(li) = integral_vec_a_r(li) - (integral2(cost_fun1_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun1_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + integral2(cost_fun1_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4)));
    if lambda > 0
        integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * integral2(cost_fun2_hdl, X_u_coord(1), X_u_coord(2), X_u_coord(3), X_u_coord(4));
    end
end
obj_sos = obj_sos + c_a_sym(1:Qa)' * integral_vec_a_r;

integral_vec_w_r = zeros(Qw, 1);
for li=1:Qw
    cost_fun2_hdl = matlabFunction(Psi_sym(li) / (poly_b_sym^Alph));
    integral_vec_w_r(li) = integral_vec_w_r(li) - integral2(cost_fun2_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) + ... 
                                                                  integral2(cost_fun2_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + ...
                                                                  integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + ...
                                                                  integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4));
%     integral_vec_w_r(li) = integral_vec_w_r(li) - (integral2(cost_fun2_hdl, X_r_coord(1), X_excld(1), X_r_coord(3), X_r_coord(4)) +  integral2(cost_fun2_hdl, X_excld(2), X_r_coord(2), X_r_coord(3), X_r_coord(4)) + integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_r_coord(3), X_excld(3)) + integral2(cost_fun2_hdl, X_r_coord(1), X_r_coord(2), X_excld(4), X_r_coord(4)));
end
obj_sos = obj_sos + c_w_sym(1:Qw)' * integral_vec_w_r;
sos_prog = sossetobj(sos_prog, obj_sos);


%%

syms x1 x2
x = [x1; x2];
% constraint_sym = poly_b_sym*(1+Alph)*sum(diag(jacobian(F*poly_a_sym+G*poly_c_sym,x))) - Alph*sum(diag(jacobian(poly_b_sym*(poly_a_sym*F+G*poly_c_sym))));
% constraint_sym

% term1
term(1) = vec_div(F, x)*poly_a_sym + jacobian(poly_a_sym, x)*F;
term(2) = vec_div(G, x)*poly_c_sym + jacobian(poly_c_sym, x) * G;
term(3) = vec_div(F, x)*poly_ab_sym + jacobian(poly_ab_sym, x) * F;
term(4) = vec_div(G, x)*poly_bc_sym + jacobian(poly_bc_sym, x) * G;

constraint_sym = (1+Alph)*poly_b_sym*sum(term(1:2)) - Alph*sum(term(3:4));

constraint_sym = vpa(constraint_sym, 8);
sos_prog = sosineq(sos_prog, constraint_sym);
sos_prog = sosineq(sos_prog, poly_a_sym);
% 
% c_s1 = mpvar('c_s1', [Qs1, 1]);
% c_s1 = [c_s1(1:Qs1); zeros(total_deg-Qs1,1)];
% s1 = p2s(c_s1' * Psi);
% c_s1_sym = p2s(c_s1);
% sos_prog = sosdecvar(sos_prog, c_s1_sym(1:Qs1));
% sos_prog = sosineq(sos_prog, s1);
% % 
% c_s3 = mpvar('c_s3', [Qs1, 1]);
% c_s3 = [c_s3(1:Qs1); zeros(total_deg-Qs1,1)];
% s3 = p2s(c_s3' * Psi);
% c_s3_sym = p2s(c_s3);
% sos_prog = sosdecvar(sos_prog, c_s3_sym(1:Qs1));
% sos_prog = sosineq(sos_prog, s3);

geo = domain_definition();
X_init = geo.X_init;
X_u = geo.X_u;
X_r = geo.X_r;
X_d = geo.X_d;
X = geo.X;

% % Constraints on X_init: div(pf+gu) > 1-eps
% ineq = constraint_sym - (1 - eps) * (poly_b_sym)^(Alph+1) - s1*X_init;
% ineq = vpa(ineq, 9);
% sos_prog = sosineq(sos_prog, ineq);

% % constraints on X_u:div(pf+gu) < eps
% ineq3 = - constraint_sym + eps*(poly_b_sym)^(Alph+1) - s3*X_u;
% ineq3 = vpa(ineq3, 9);
% sos_prog = sosineq(sos_prog, ineq3);


% c_s2 = mpvar('c_s2', [Qs1, 1]);
% c_s2 = [c_s2(1:Qs1); zeros(total_deg-Qs1,1)];
% s2 = p2s(c_s2' * Psi);
% c_s2_sym = p2s(c_s2);
% sos_prog = sosdecvar(sos_prog, c_s2_sym(1:Qs1));
% 
% sos_prog = sosineq(sos_prog, s2);
% ineq2 = - poly_a_sym + (-1e-5*(poly_b_sym)^(Alph+1)) - s2*X_u;
% ineq2 = vpa(ineq2, 9);
% sos_prog = sosineq(sos_prog, ineq2);

% Constraints on X_init: div(pf+gu) > 1-eps on X_r
% ineq = poly_a_sym - (1-1e-5)*(poly_b_sym)^(Alph+1) - s1*X_r;
% ineq = vpa(ineq, 9);
% sos_prog = sosineq(sos_prog, ineq);

%% solve
disp('Solving sos program')

option.solver = 'sedumi';
option.params.tol = 1e-15;
[sos_prog, info] = sossolve(sos_prog, option);

disp('==================Solution check for the optimal control problem================')
disp('----------- feasibility check --------------')
disp('info.pinf')
info.pinf
disp('info.dinf')
info.dinf
if info.pinf == 1
    disp('Primal infeasibility')
elseif info.dinf == 1
    disp('Dual infeasibility')
else
    disp('------------- Found feasible safety verification rho function ---------------')
end

disp('------- numerical issue -----------')
info.numerr
if info.numerr == 2
    disp('Attention! numerical issues! Do not trust the SDP results')
    keyboard;
elseif info.numerr == 1
    disp('Attention! numerical inaccuracy!')
end

%%
% get solution: cx, ax. 
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym);
cx_sol = s2p(cx_sol_sym);
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym);
ax_sol = s2p(ax_sol_sym);
wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
wx_sol = s2p(wx_sol_sym);
bx_sol_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol = s2p(bx_sol_sym);

ax_sol_sym = vpa(ax_sol_sym, 8)
bx_sol_sym = vpa(bx_sol_sym, 8);
cx_sol_sym = vpa(cx_sol_sym, 8)
wx_sol_sym = vpa(wx_sol_sym, 8);

%% safety verification

% constraints for safety
% 1. rho > 0 in X_init
ineq_safe_1 = poly_safe_rho_sym - poly_safe_s1 * X_init;
ineq_safe_1 = vpa(ineq_safe_1, 8);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_1);

% 2. div(rho*f + gu) >= 0 in X - X_r
[F, G] = dynamics_vdp();
F = p2s(F);

div_pf_sym = vec_div(poly_safe_rho_sym.*F, [x1, x2]);
num_div_gu = diff(poly_safe_rho_sym*cx_sol_sym, x2)*ax_sol_sym - transpose((G.*cx_sol_sym*poly_safe_rho_sym)) * vec_grad(ax_sol_sym, [x1, x2]);
num_div_gu = vpa(num_div_gu, 8);
div_pf_sym = vpa(div_pf_sym, 8);
ineq_safe_2 = div_pf_sym * ax_sol_sym^2 + num_div_gu - poly_safe_s2 * X + poly_safe_s3 * X_r  ;
ineq_safe_2 = vpa(ineq_safe_2, 8);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_2);

% 3. rho < 0 in X_u + (X_d - X)
ineq_safe_3 = - poly_safe_rho_sym - poly_safe_s6*X_u;
ineq_safe_3 = vpa(ineq_safe_3, 8);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_3);

ineq_safe_4 = - poly_safe_rho_sym - poly_safe_s4*X_d + poly_safe_s5*X;
ineq_safe_4 = vpa(ineq_safe_4, 8);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_4);

% 4. rho > 0 in X_r
ineq_safe_5 = poly_safe_rho_sym - poly_safe_s7 * X_r;
ineq_safe_5 = vpa(ineq_safe_5, 8);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_5);


%% solve the safety verifycation feasibility problem
option.solver = 'sedumi';
option.params.tol = 1e-15;
[sos_safety_prog, info] = sossolve(sos_safety_prog, option);

%% solutions
disp('============Solution check for the safety verification problem=============')
disp('----------- feasibility check --------------')
disp('info.pinf')
info.pinf
disp('info.dinf')
info.dinf
if info.pinf == 1
    disp('Primal infeasibility')
elseif info.dinf == 1
    disp('Dual infeasibility')
else
    disp('------------- Found feasible safety verification rho function ---------------')
end

disp('------- numerical issue -----------')
disp('info.numerr')
info.numerr
if info.numerr == 2
    disp('Attention! numerical issues! Do not trust the SDP results')
    keyboard;
elseif info.numerr == 1
    disp('Numerical inaccuracy!')
end

disp('------------ Solution for Rho---------------')
rho_safe_sol = sosgetsol(sos_safety_prog, poly_safe_rho_sym)

if rho_safe_sol == 0
    disp('Trivial solution of safety-reachable Rho')
end

disp('Saving solved data: ')
name = ['experiments/results_SOS_Rational_lbd_',num2str(lambda),'_eps_',num2str(eps),'.mat'];
disp(name)
save(name)

