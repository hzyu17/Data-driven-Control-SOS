function [data_file_name] = solve_optimal_control(lambda, eps, domain_geo)
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
ax_sol_sym = p2s(clean_polynomial(s2p(ax_sol_sym), trunc_decim))
ax_sol_sym = vpa(ax_sol_sym, truncation_precision);
ax_sol = s2p(ax_sol_sym);

wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
wx_sol_sym = vpa(wx_sol_sym, truncation_precision);
wx_sol_sym = p2s(clean_polynomial(s2p(wx_sol_sym), trunc_decim));
wx_sol_sym = vpa(wx_sol_sym, truncation_precision)
wx_sol = s2p(wx_sol_sym);

bx_sol_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol_sym = p2s(clean_polynomial(s2p(bx_sol_sym), trunc_decim));
bx_sol_sym = vpa(bx_sol_sym, truncation_precision)
bx_sol = s2p(bx_sol_sym);

for i_c = 1:dim_m
    cx_sol_sym(i_c) = vpa(cx_sol_sym(i_c), truncation_precision);
end
 
 %% saving data
disp('----- Saving solved optimal control data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_file_name = ['experiments/',date,'optimal_control_results_SOS_Rational_lbd_',num2str(lambda),'.mat'];
disp(data_file_name)
save(data_file_name)

