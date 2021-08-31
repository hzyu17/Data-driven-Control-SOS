function [data_file_name] = ocp_safety_rational_main(ocp_data_file)
%% SOS program creation or reading
disp('Loading predefined sos program')
load('sos_prog.mat')

%% calculate the OCP objective function (linear function of c_a and c_c)
% Psi_sym = p2s(Psi);
num_xu = domain_geo.num_xu;

integral_vec_a_r = zeros(Qa, 1);

% TODO: modulate this integration paragraph, and write code for multiple
% unsafe sets and initial sets. Done for the multiple unsafe set case.
% 06/30: moved the whole integration part in the program construction part.
% Done with the multiple unsafe sets case.
if lambda > 0
    for li=1:Qa
        for i_xu = 1: num_xu % multiple unsafe sets
            X_u_coord_i = domain_geo.X_u_coord(i_xu, :);
            integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * integral2(cost_fun2_hdl, X_u_coord_i(1), X_u_coord_i(2), ...
                                                                                                                                                                                 X_u_coord_i(3), X_u_coord_i(4));
        end
    end
end
% transpose(c_a_sym(1:Qa)) * integral_vec_a_r
obj_sos = obj_sos + transpose(c_a_sym(1:Qa)) * integral_vec_a_r;

%% setting the objective for the optimal control problem
obj_sos
obj_sos = clean_polynomial(obj_sos, 1e-9);
obj_sos = vpa(obj_sos, 8);
obj_sos
sos_prog = sossetobj(sos_prog, obj_sos);
% sos_prog = sossetobj(sos_prog, 1);

%% program definition
syms x1 x2
x = [x1; x2];
truncation_precision = 9;
% truncation_precision = 9
% trunc_decim = 1e-12;
trunc_decim = 1e-9;

% constraint_sym = poly_b_sym*(1+Alph)*sum(diag(jacobian(F*poly_a_sym+G*poly_c_sym,x))) - Alph*sum(diag(jacobian(poly_b_sym*(poly_a_sym*F+G*poly_c_sym))));
% constraint_sym

% terms
% term(1) = vec_div(F, x)*poly_a_sym + jacobian(poly_a_sym, x)*F;
% term(2) = vec_div(G, x)*poly_c_sym + reshape(G, 1, []) * reshape(transpose(jacobian(poly_c_sym, x)), [], 1);
% term(3) = vec_div(F, x)*poly_ab_sym + jacobian(poly_ab_sym, x) * F;
% term(4) = vec_div(G, x)*poly_bc_sym +  reshape(G, 1, []) * reshape(transpose(jacobian(poly_bc_sym, x)), [], 1);

% multi-dimensional inputs
term(1) = vec_div(F.*poly_a_sym, x);
term(2) = vec_div(G*poly_c_sym, x);
term(3) = vec_div(F.*poly_ab_sym, x);
term(4) = vec_div(G*poly_bc_sym, x);

constraint_sym = (1+Alph)*poly_b_sym*sum(term(1:2)) - Alph*sum(term(3:4));

constraint_sym = vpa(constraint_sym, truncation_precision);
sos_prog = sosineq(sos_prog, constraint_sym);
sos_prog = sosineq(sos_prog, poly_a_sym);

%% solve
disp('------------- Solving optimal control sos program -------------')

option.solver = 'sedumi';
option.params.tol = 1e-15; % tolerance for the stop
% option.params.bigeps = 1e-9; % the tolerence for the dual infeasibility check 
% option.params.vplot = 0;
% option.params.free = 0;
% option.params.stepdif = 1;
% option.params.alg = 1;
option.params.errors = 1;

sos_prog = solve_show_info(sos_prog, option);

%%
% get solution: cx, ax.
disp('----- get solution from the sos program -----')
cx_sol = [];
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym)
for i_c = 1:dim_m
    cx_sol_sym(i_c) = vpa(cx_sol_sym(i_c), truncation_precision);
    if cx_sol_sym(i_c) == 0
        keyboard
    end
    cx_sol = [cx_sol; s2p(cx_sol_sym(i_c))];
end
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym);
ax_sol_sym = clean_polynomial(ax_sol_sym, trunc_decim)
ax_sol_sym = vpa(ax_sol_sym, truncation_precision);
ax_sol = s2p(ax_sol_sym);

wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
wx_sol_sym = vpa(wx_sol_sym, truncation_precision);
wx_sol_sym = clean_polynomial(wx_sol_sym, trunc_decim);
wx_sol_sym = vpa(wx_sol_sym, truncation_precision)
wx_sol = s2p(wx_sol_sym);

bx_sol_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol_sym = clean_polynomial(bx_sol_sym, trunc_decim);
bx_sol_sym = vpa(bx_sol_sym, truncation_precision)
bx_sol = s2p(bx_sol_sym);

for i_c = 1:dim_m
    cx_sol_sym(i_c) = vpa(cx_sol_sym(i_c), truncation_precision);
end

%% safety verification
disp('----- the safety eventuality verification problem -----')
% constraints for safety
% 1. rho > 0 in X_init
ineq_safe_1 = poly_safe_rho_sym - poly_safe_s1 * domain_geo.poly_X_init;
% ineq_safe_1 = p2s(clean_polynomial(s2p(ineq_safe_1), trunc_decim));
% ineq_safe_1 = vpa(ineq_safe_1, truncation_precision);
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_1);

% 2. div(rho*(f + gu)) >= 0 in X - X_r
[F, G] = dynamics_definition();

if strcmp(class(F), 'polynomial')
    F = p2s(F);
end

div_pf_sym = vec_div(poly_safe_rho_sym.*F, [x1, x2]);
num_div_gu = vec_div(G*cx_sol_sym.*poly_safe_rho_sym, [x1;x2])*ax_sol_sym - ...
                               transpose((G*cx_sol_sym.*poly_safe_rho_sym)) * vec_grad(ax_sol_sym, [x1, x2]);
% num_div_gu = p2s(clean_polynomial(s2p(num_div_gu), trunc_decim));

ineq_safe_2 = div_pf_sym * ax_sol_sym^2 + num_div_gu + poly_safe_s3 * domain_geo.poly_X_r - poly_safe_s2 * domain_geo.poly_X;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_2);

% 3. rho <= 0 in X_u + (X_d - X)
ineq_safe_3 = - poly_safe_rho_sym - poly_safe_s6*domain_geo.poly_X_u;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_3);

ineq_safe_4 = - poly_safe_rho_sym - poly_safe_s4*domain_geo.poly_X_d + poly_safe_s5*domain_geo.poly_X;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_4);

%% solve the safety verifycation feasibility problem
disp('----- solving the safety eventuality problem -----')
option_safety.solver = 'sedumi';
option_safety.params.tol = 1e-15;
% option.params.free = 0;
% option.params.stepdif = 1;
% option.params.vplot = 1;

% solve and show info
sos_safety_prog = solve_show_info(sos_safety_prog, option_safety);

disp('------------ Solution for Rho---------------')
rho_safe_sol = sosgetsol(sos_safety_prog, poly_safe_rho_sym);
rho_safe_sol = clean_polynomial(rho_safe_sol, trunc_decim)
rho_safe_sol = vpa(rho_safe_sol, truncation_precision);

if rho_safe_sol == 0
    disp('Trivial solution of safety-reachable Rho')
end

%% manually verify the solved polynomials
% optimal control
ab_sol_sym = ax_sol_sym * poly_b_sym;
bc_sol_sym = cx_sol_sym * poly_b_sym;

sol_term(1) = vec_div(F.*ax_sol_sym, x);
sol_term(2) = vec_div(G*cx_sol_sym, x);
sol_term(3) = vec_div(F.*ab_sol_sym, x);
sol_term(4) = vec_div(G*bc_sol_sym, x);

constraint_sym_sol = (1+Alph)*poly_b_sym*sum(sol_term(1:2)) - Alph*sum(sol_term(3:4));

 res1 = findsos(constraint_sym_sol);
 res2 = findsos(ax_sol_sym);
 
%  safety program
 safe_s1_sol = sosgetsol(sos_safety_prog, poly_safe_s1);
 safe_s2_sol = sosgetsol(sos_safety_prog, poly_safe_s2);
 safe_s3_sol = sosgetsol(sos_safety_prog, poly_safe_s3);
 safe_s4_sol = sosgetsol(sos_safety_prog, poly_safe_s4);
 safe_s5_sol = sosgetsol(sos_safety_prog, poly_safe_s5);
 safe_s6_sol = sosgetsol(sos_safety_prog, poly_safe_s6);
 
 ineq_safe_1_sol = sosgetsol(sos_safety_prog, ineq_safe_1);
 ineq_safe_2_sol = sosgetsol(sos_safety_prog, ineq_safe_2);
 ineq_safe_3_sol = sosgetsol(sos_safety_prog, ineq_safe_3);
 ineq_safe_4_sol = sosgetsol(sos_safety_prog, ineq_safe_4);

%  ineq_safe_1_sol = rho_safe_sol - safe_s1_sol * domain_geo.poly_X_init;
%  
%  div_pf_sym = vec_div(rho_safe_sol.*F, [x1, x2]);
%  num_div_gu_sol = vec_div(G*cx_sol_sym.*rho_safe_sol, [x1;x2])*ax_sol_sym - transpose((G*cx_sol_sym.*rho_safe_sol)) * ...
%      vec_grad(ax_sol_sym, [x1, x2]);
%  num_div_gu_sol = vpa(num_div_gu_sol, 8);
%  div_pf_sym = vpa(div_pf_sym, 8);
%  ineq_safe_2_sol = div_pf_sym * ax_sol_sym^2 + num_div_gu_sol - safe_s2_sol * domain_geo.poly_X + safe_s3_sol * domain_geo.poly_X_r;
%  
%  ineq_safe_3_sol = - rho_safe_sol - safe_s6_sol*domain_geo.poly_X_u;
%  
%  ineq_safe_4_sol = - rho_safe_sol - safe_s4_sol*domain_geo.poly_X_d + safe_s5_sol*domain_geo.poly_X;
 
 res3 = findsos(safe_s1_sol);
 res4 = findsos(safe_s2_sol);
 res5 = findsos(safe_s3_sol);
 res6 = findsos(safe_s4_sol);
 
 res7 = findsos(ineq_safe_1_sol);
 res8 = findsos(ineq_safe_2_sol);
 res9 = findsos(ineq_safe_3_sol);
 res10 = findsos(ineq_safe_4_sol);
 
 if isempty(res1) || isempty(res2) || isempty(res3) || isempty(res4) || isempty(res5) || isempty(res6) || isempty(res7) || isempty(res8) || isempty(res9) || isempty(res10)
 disp('Infeasible solution!')
 keyboard
 end
 
 %% saving data
disp('----- Saving solved data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_file_name = ['experiments/',date,'_results_SOS_Rational_lbd_',num2str(lambda),'.mat'];
disp(data_file_name)
save(data_file_name)

