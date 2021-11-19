function [data_file_name] = solve_optimal_control(lambda, sos_prog_file, gEDMD_file, option)
%% preparations
disp('------ Loading predefined sos program -------')

% load file
load(sos_prog_file);

% data driven or not
data_approx_option = option.data_driven_option;

% geometry
geometry = domain_definition();

%numerical truncations
solution_truncation = option.solution_truncation;
trunc_decim = option.variable_trunc;

%% objective function (linear function of c_a and c_c)
% Psi_sym = p2s(Psi);
num_xu = geometry.num_xu;
integral_vec_a_r = zeros(Qa, 1);
% TODO: modulate this integration paragraph, and write code for multiple
% unsafe sets and initial sets. Done for the multiple unsafe set case.
% 06/30: moved the whole integration part in the program construction part.
% Done with the multiple unsafe sets case.
base_integral_obj = obj_sos;
if lambda > 0
    for i_xu = 1: num_xu % multiple unsafe sets
        X_u_coord_i = geometry.X_u_coord(i_xu);
        lower_bd_Xu = matlabFunction(geometry.X_u_cen(i_xu, 2) - sqrt(geometry.r_xu(i_xu)^2 - (var_sym.x(1) - geometry.X_u_cen(i_xu, 1)).^2));
        upper_bd_Xu = matlabFunction(geometry.X_u_cen(i_xu, 2) + sqrt(geometry.r_xu(i_xu)^2 - (var_sym.x(1) - geometry.X_u_cen(i_xu, 1)).^2));
         for li=1:Qa
                 cost_fun2_hdl = matlabFunction(Psi_sym(li) ./ (poly_b_sym.^Alph));
                 % cube
                 integral_vec_a_r(li) = integral_vec_a_r(li) + lambda * integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, ... 
                                                                                                                                                      X_u_coord_i.ymin, X_u_coord_i.ymax);
                 % circle
%                 integral_vec_a_r(li) = lambda * integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, lower_bd_Xu, upper_bd_Xu);
%                 integral_vec_a_r(li) = integral_vec_a_r(li) + integral2(cost_fun2_hdl, X_u_coord_i.xmin, X_u_coord_i.xmax, X_u_coord_i.ymin, X_u_coord_i.ymax);
         end
    end
    obj_sos = obj_sos + transpose(c_a_sym(1:Qa)) * integral_vec_a_r;
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
if strcmp(data_approx_option.type, 'data_driven')
    if strcmp(data_approx_option.approx, 'gEDMD')
    %     load('gEDMD_res.mat')
%         load('sampling_data/gEDMD_res_01_15_21_01_27_58.mat') %bx = LQR
        load(gEDMD_file)
    else
        load('sampling_data/sampling_data.mat')
        load('sampling_data/EDMD_res.mat')
    end
    
    %     pvar x1 x2
    %     x = [x1;x2];
    Div_F = 0; Div_G = 0;
    F = c_x'*L{1}*Psi;
    G = [];
    
     for i1 = 1 : nx
            Div_Fi = diff(F(i1), var_poly.x(i1));
            Div_F = Div_F + Div_Fi;
     end
    Div_F.coefficient(find(abs(Div_F.coefficient) <= trunc_decim)) = 0;
    % term1
    term(1) = c_a'*L{1}*Psi + c_a'*Psi*Div_F;
    % term3
    term(3) = c_ab'*L{1}*Psi + c_ab'*Psi*Div_F;
    
    for i_input = 2:dim_m+1
        if strcmp(data_approx_option.approx, 'gEDMD')
            G = c_x'*(L{i_input}-L{1})*Psi;
        else %EDMD
            G = c_x'*L{i_input}*Psi;
        end

        for i1 = 1 : nx
            Div_Gi = diff(G(i1), var_poly.x(i1));
            Div_G = Div_G + Div_Gi;
        end
        Div_G.coefficient(find(abs(Div_G.coefficient) <= trunc_decim)) = 0;

        % term2, term4
        if strcmp(data_approx_option.approx, 'gEDMD')
%             if length(c_c) == 1
%                 c_c = c_c{1};
%             end
            if i_input == 2
                term(2) = c_c{i_input-1}'*(L{i_input}-L{1})*Psi + c_c{i_input-1}'*Psi*Div_G;
                term(4) = c_bc{i_input-1}'*(L{i_input}-L{1})*Psi + c_bc{i_input-1}'*Psi*Div_G;
            else
                term(2) = term(2) + c_c{i_input-1}'*(L{i_input}-L{1})*Psi + c_c{i_input-1}'*Psi*Div_G;
                term(4) = term(4) + c_bc{i_input-1}'*(L{i_input}-L{1})*Psi + c_bc{i_input-1}'*Psi*Div_G;
            end
        else %EDMD
            if i_input == 2
                term(2) = c_c{i_input-1}'*L{i_input}*Psi + c_c{i_input-1}'*Psi*Div_G;
                term(4) = c_bc{i_input-1}'*L{i_input}*Psi + c_bc{i_input-1}'*Psi*Div_G;
            else
                term(2) = term(2) + c_c{i_input-1}'*L{i_input}*Psi + c_c{i_input-1}'*Psi*Div_G;
                term(4) = term(4) + c_bc{i_input-1}'*L{i_input}*Psi + c_bc{i_input-1}'*Psi*Div_G;
            end
        end
    end
    constraint_poly = (1+Alph).*poly_b*sum(term(1:2)) - Alph.*sum(term(3:4));
    % L1-DATA-DRIVEN
    %     constraint_poly.coefficient(find(abs(constraint_poly.coefficient)<=1e-1)) = 0;
    constraint_poly.coefficient(find(abs(constraint_poly.coefficient)<=solution_truncation)) = 0;
    constraint_sym = p2s(constraint_poly);
elseif strcmp(data_approx_option.type, 'model_based')
    term(1) = vec_div(dynamics_sym.F.*poly_a_sym, var_sym.x);
    term(2) = vec_div(dynamics_sym.G*poly_c_sym, var_sym.x);
    term(3) = vec_div(dynamics_sym.F.*poly_ab_sym, var_sym.x);
    term(4) = vec_div(dynamics_sym.G*poly_bc_sym, var_sym.x);
    constraint_sym = (1+Alph).*poly_b_sym.*sum(term(1:2)) - Alph.*sum(term(3:4));
end

% vdp
% constraint_sym = clean_polynomial(constraint_sym, trunc_decim);
% constraint_sym = vpa(constraint_sym, truncation_precision);

sos_prog = sosineq(sos_prog, constraint_sym);

%% constraint on initial set X_0
geo = domain_definition();
X_0 = geo.poly_X_init;
alhpa_X0 = option.optimization.alpha_x0;
constaint_x0_i = constraint_sym-alhpa_X0*X_0(1);
sos_prog = sosineq(sos_prog, constaint_x0_i);

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

%% operations on the solution
 options_solutions.solver = 'sedumi';
 
% get solution: cx, ax.
disp('----- get solution from the sos program -----')
cx_sol = [];
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym);
for i_c = 1:dim_m
%     cx_sol_sym(i_c) = clean_polynomial(cx_sol_sym(i_c), trunc_decim);
    if cx_sol_sym(i_c) == 0
        keyboard
    end
    cx_sol_sym_i = cx_sol_sym(i_c);
    cx_sol = [cx_sol; s2p(cx_sol_sym(i_c))];
end
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym);
% coeff_ax = get_coefficients(ax_sol_sym);
% vdp
% ax_sol_sym = clean_polynomial(ax_sol_sym, trunc_decim);
% ax_sol_sym = vpa(ax_sol_sym, truncation_precision);
ax_sol = s2p(ax_sol_sym);

% find sos for a(x)
disp('----- Is the solved w(x) positive definite? -----')
[Q_sol_a, ~] = findsos(ax_sol_sym, [], options_solutions);
eig_vals = eig((Q_sol_a + Q_sol_a')/2); 
is_psd_sol_a = [sum([sign(eig_vals)==1])==length(Q_sol_a)]
if ~is_psd_sol_a
    keyboard
end

% poly_b_sym = sosgetsol(sos_prog, poly_b_sym);
bx_sol = s2p(poly_b_sym);

%% L2 case
% check the psd for the matrix M and w(x)
if strcmp(option.control_penalty_type, 'L2')
    % find sos for M
    M_sol = sosgetsol(sos_prog,  M);    
    disp('----- Is the solved M positive definite? -----')
  
    [Q_sol_M, ~, ~] = findsos(M_sol, [], options_solutions);
    eig_vals = eig((Q_sol_M + Q_sol_M')/2); 
    is_pd_sol_M = [sum([sign(eig_vals)==1])==length(Q_sol_M)]
    if ~is_pd_sol_M
        keyboard
    end
    
    % find sos for w(x)
    wx_sol_sym = sosgetsol(sos_prog, poly_w_sym);
    wx_sol = s2p(wx_sol_sym);
    
    disp('----- Is the solved w(x) positive definite? -----')
    [Q_sol_w, ~] = findsos(wx_sol_sym, [], options_solutions);
    eig_vals = eig((Q_sol_w + Q_sol_w')/2); 
    is_pd_sol_w = [sum([sign(eig_vals)==1])==length(Q_sol_w)]
    if ~is_pd_sol_w
        keyboard
    end
  
    % evaluate the solved objective function
    %% Base integral
    % monomial integrations
    geometry = domain_definition();
    lower_bd_X = matlabFunction(geometry.X_cen(2) - sqrt(geometry.r_X^2 - (var_sym.x(1) - geometry.X_cen(1)).^2));
    upper_bd_X = matlabFunction(geometry.X_cen(2) + sqrt(geometry.r_X^2 - (var_sym.x(1) - geometry.X_cen(1)).^2));

    lower_bd_Xr = matlabFunction(geometry.X_r_cen(2) - sqrt(geometry.r_xr^2 - (var_sym.x(1) - geometry.X_r_cen(1)).^2));
    upper_bd_Xr = matlabFunction(geometry.X_r_cen(2) + sqrt(geometry.r_xr^2 - (var_sym.x(1) - geometry.X_r_cen(1)).^2));


    eval_cost_fun_hdl = matlabFunction(ax_sol_sym.* poly_q_sym ./ (poly_b_sym.^Alph) + transpose(cx_sol_sym)* R_ocp * cx_sol_sym /  ax_sol_sym);
    eval_cost_fun_hdl1 = matlabFunction(ax_sol_sym.* poly_q_sym ./ (poly_b_sym.^Alph) + wx_sol_sym);
% circle
    disp('L2 cost function recalculated')
    eval_cost_fun = integral2(eval_cost_fun_hdl, geometry.X_0.xmin, geometry.X_r_coord.xmin, lower_bd_X, upper_bd_X) + ...
                                             integral2(eval_cost_fun_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, lower_bd_X, lower_bd_Xr) + ...
                                             integral2(eval_cost_fun_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, upper_bd_Xr, upper_bd_X) + ...
                                             integral2(eval_cost_fun_hdl, geometry.X_r_coord.xmax, geometry.X_0.xmax, lower_bd_X, upper_bd_X)
                                         
    disp('L2 cost function recalculated using wx')
    eval_cost_fun = integral2(eval_cost_fun_hdl1, geometry.X_0.xmin, geometry.X_r_coord.xmin, lower_bd_X, upper_bd_X) + ...
                                             integral2(eval_cost_fun_hdl1, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, lower_bd_X, lower_bd_Xr) + ...
                                             integral2(eval_cost_fun_hdl1, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, upper_bd_Xr, upper_bd_X) + ...
                                             integral2(eval_cost_fun_hdl1, geometry.X_r_coord.xmax, geometry.X_0.xmax, lower_bd_X, upper_bd_X)
                                         
end

% evaluate the cost function
disp('obj_sos_sol')
obj_sos_sol = sosgetsol(sos_prog, obj_sos)
disp('base_obj_sos_sol')
base_obj_sos_sol = sosgetsol(sos_prog,  base_integral_obj)
constraint_sym_sol = sosgetsol(sos_prog, constraint_sym)

date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_driven_type = option.data_driven_option.type;
regularizing_type = option.control_penalty_type;
file_prefix = ['res/',date, '_', dynamics_option,'_trj_',regularizing_type,'_',data_driven_type,'_','lamb_',num2str(lambda),'ocp_feas_',num2str(ocp_info.feasratio)];
name_div_surface = [file_prefix,'_constraintsym_surface_alpha_',num2str(alhpa_X0),'.png'];
name_h0_surface = [file_prefix,'_h0_surface_alpha_',num2str(alhpa_X0),'.png'];

fig1 = figure();
plot_sets(geo);
fsurf(constraint_sym_sol, 'FaceAlpha', 0.3);
xlim([-4.5, 4.5]);
ylim([-4.5, 4.5]);
view(-19,56);

fig2 = figure();
plot_sets(geo);
fsurf(constraint_sym_sol/poly_b_sym, 'FaceAlpha', 0.3);
xlim([-4.5, 4.5]);
ylim([-4.5, 4.5]);
view(-19,56);

saveas(fig1, name_div_surface);
saveas(fig2, name_h0_surface);

 %% data file
 data_driven_type = option.data_driven_option.type;
 input_regularizer = option.control_penalty_type;
disp('----- Saving solved optimal control data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));

data_file_name = ['experiments/',date,'_',dynamics_option,...
    '_optimal_control_lbd_', num2str(lambda), '_',input_regularizer,'_',data_driven_type,'_feasratio_', num2str(ocp_info.feasratio),'.mat'];
% disp(data_file_name);

% save file
save(data_file_name);

%% plotting elements
% save in the mat file
plotting_elements = plotting_handler(data_file_name, option);

%% save file:  added the plotting elements
save(data_file_name)

