function [data_file_name] = solve_optimal_control(lambda, sos_prog_file)
%% preparations
disp('------ Loading predefined sos program -------')
sos_prog_file
load(sos_prog_file);

geometry = domain_definition();

% numerical truncation
% truncation_precision = 12;
truncation_precision = 9;
% trunc_decim = 1e-12;
trunc_decim = 1e-12;

%% objective function (linear function of c_a and c_c)
% Psi_sym = p2s(Psi);
num_xu = geometry.num_xu;
integral_vec_a_r = zeros(Qa, 1);
% TODO: modulate this integration paragraph, and write code for multiple
% unsafe sets and initial sets. Done for the multiple unsafe set case.
% 06/30: moved the whole integration part in the program construction part.
% Done with the multiple unsafe sets case.

if lambda > 0
    for i_xu = 1: num_xu % multiple unsafe sets
        X_u_coord_i = geometry.X_u_coord(i_xu);
        lower_bd_Xu = matlabFunction(geometry.X_u_cen(i_xu, 2) - sqrt(geometry.r_xu(i_xu)^2 - (var_sym.x(1) - geometry.X_u_cen(i_xu, 1)).^2));
        upper_bd_Xu = matlabFunction(geometry.X_u_cen(i_xu, 2) + sqrt(geometry.r_xu(i_xu)^2 - (var_sym.x(1) - geometry.X_u_cen(i_xu, 1)).^2));
         for li=1:Qa
                 cost_fun2_hdl = matlabFunction(Psi_sym(li) ./ (poly_b_sym.^Alph));
                 integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, ... 
                                                                                                                                                      X_u_coord_i.ymin, X_u_coord_i.ymax);
%                 integral_vec_a_r(li) = lambda * integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, lower_bd_Xu, upper_bd_Xu)
%                 integral_vec_a_r(li) = integral_vec_a_r(li) + integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, X_u_coord_i.ymin, X_u_coord_i.ymax);
         end
    end
    obj_sos = obj_sos + transpose(c_a_sym(1:Qa)) * integral_vec_a_r
end

 % rescale objective function to provide a better numerical condition
%  obj_sos = obj_sos / (1e2*min(abs(get_coefficients(obj_sos))));

%% setting the objective for the optimal control problem
sos_prog = sossetobj(sos_prog, obj_sos);

%% program definition
% multi-dimensional inputs
if exist('dynamics_sym')==0
    dynamics_sym = dynamics_definition('sym', dynamics_option);
end

term(1) = vec_div(dynamics_sym.F.*poly_a_sym, var_sym.x);
term(2) = vec_div(dynamics_sym.G*poly_c_sym, var_sym.x);
term(3) = vec_div(dynamics_sym.F.*poly_ab_sym, var_sym.x);
term(4) = vec_div(dynamics_sym.G*poly_bc_sym, var_sym.x);

constraint_sym = (1+Alph).*poly_b_sym.*sum(term(1:2)) - Alph.*sum(term(3:4));

% vdp
% constraint_sym = clean_polynomial(constraint_sym, trunc_decim);
% constraint_sym = vpa(constraint_sym, truncation_precision);

sos_prog = sosineq(sos_prog, constraint_sym);
sos_prog = sosineq(sos_prog, poly_a_sym);

%% solve
disp('------------- Solving optimal control sos program -------------')
option.solver = 'sedumi';
option.params.tol = 1e-15; % tolerance for the stop
% option.params.stepdif = 1;
% option.params.vplot = 1;
% option.params.alg = 1;
% option.params.theta = 1;
% option.params.beta = 0.9;

[sos_prog, ocp_info] = solve_show_info(sos_prog, option);

%%
% get solution: cx, ax.
disp('----- get solution from the sos program -----')
cx_sol = [];
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym);
for i_c = 1:dim_m
%     cx_sol_sym(i_c) = clean_polynomial(cx_sol_sym(i_c), trunc_decim);
    if cx_sol_sym(i_c) == 0
        keyboard
    end
    cx_sol_sym_i = cx_sol_sym(i_c)
    cx_sol = [cx_sol; s2p(cx_sol_sym(i_c))];
end
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym);
% coeff_ax = get_coefficients(ax_sol_sym);

% vdp
% ax_sol_sym = clean_polynomial(ax_sol_sym, trunc_decim);
% ax_sol_sym = vpa(ax_sol_sym, truncation_precision);
ax_sol = s2p(ax_sol_sym)

wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
% coeff_wx = get_coefficients(wx_sol_sym);
% wx_sol_sym = p2s(clean_polynomial(s2p(wx_sol_sym), trunc_decim));
% wx_sol_sym = vpa(wx_sol_sym, truncation_precision)
wx_sol = s2p(wx_sol_sym);

bx_sol_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol = s2p(bx_sol_sym);
 
 %% data file
disp('----- Saving solved optimal control data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_file_name = ['experiments/',date,'_',dynamics_option,...
    '_optimal_control_lbd_', num2str(lambda), '_',input_regularizer,'_feasratio_', num2str(ocp_info.feasratio),'.mat'];
% disp(data_file_name);

% save file
save(data_file_name);

%% plotting elements
% save in the mat file
plotting_elements = plotting_handler(data_file_name);

%% save file:  added the plotting elements
save(data_file_name)

