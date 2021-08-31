function data_file_name = solve_safety_verification(ocp_data_file, safety_prog_file)
%% safety verification
load(ocp_data_file, 'cx_sol_sym', 'ax_sol_sym', 'nx', 'x1', 'x2', 'dynamics_sym');

%% create safety program
load('polynomials_def.mat')
geometry = domain_definition();

% define sos program
sos_safety_prog = sosprogram(var_sym.x);

% degrees
% VDP
% deg_safe_rho = 11;
deg_safe_rho = 12;

% deg_safe_rho = 3;
Q_safe_rho = nchoosek(deg_safe_rho+nx, nx);
c_safe_rho = mpvar('c_safe_rho', [Q_safe_rho, 1]);
c_safe_rho = [c_safe_rho(1:Q_safe_rho); zeros(nPsi-Q_safe_rho, 1)];
c_safe_rho_sym = p2s(c_safe_rho);
poly_safe_rho_sym = p2s(transpose(c_safe_rho) * Psi);

sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_rho_sym(1:Q_safe_rho));
sparse = 'false';

deg_pX = polynomialDegree(geometry.poly_X);
deg_pXinit = polynomialDegree(geometry.poly_X_init);
deg_pXr = polynomialDegree(geometry.poly_X_r);
deg_pXd = polynomialDegree(geometry.poly_X_d);
deg_pXu = polynomialDegree(geometry.poly_X_u);

deg_safe_s1 = deg_safe_rho - deg_pXinit;
[sos_safety_prog, poly_safe_s1] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s1', Psi, deg_safe_s1,sparse);

% 2-dimensional case
deg_safe_s4 = deg_safe_rho - deg_pXd;
[sos_safety_prog, poly_safe_s4] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s4', Psi, deg_safe_s4,sparse);

deg_safe_s5 = deg_safe_rho - deg_pX;
[sos_safety_prog, poly_safe_s5] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s5', Psi, deg_safe_s5,sparse);

deg_safe_s6 = deg_safe_rho - deg_pXu;
[sos_safety_prog, poly_safe_s6] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s6', Psi, deg_safe_s6,sparse);

%% solve safety program
% safety_prog_file
disp('----- the safety eventuality verification problem -----')
geometry = domain_definition();

% constraints for safety
% 1. rho > 0 in X_init
% eps = 0.1;
ineq_safe_1= poly_safe_rho_sym - poly_safe_s1 * geometry.poly_X_init;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_1);

% 2. div(rho*(f + gu)) >= 0 in X - X_r
% dynamics = dynamics_definition();
F = dynamics_sym.F;
G = dynamics_sym.G;

div_pf_sym = vec_div(poly_safe_rho_sym.*F, [x1, x2]);
num_div_gu = vec_div(G*cx_sol_sym.*poly_safe_rho_sym, [x1;x2])*ax_sol_sym - ...
                               transpose((G*cx_sol_sym.*poly_safe_rho_sym)) * vec_grad(ax_sol_sym, [x1, x2]);
num_div_rhof_gu = div_pf_sym * ax_sol_sym^2 + num_div_gu;

deg_div_rhof_gu = polynomialDegree(num_div_rhof_gu);

deg_safe_s2 = deg_div_rhof_gu - deg_pX;
deg_diff = total_deg - deg_safe_s2
[sos_safety_prog, poly_safe_s2] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s2', Psi, deg_safe_s2,sparse);

deg_safe_s3 =  deg_div_rhof_gu - deg_pXr;
deg_diff = total_deg - deg_safe_s3
[sos_safety_prog, poly_safe_s3] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s3', Psi, deg_safe_s3,sparse);

% num_div_gu = p2s(clean_polynomial(s2p(num_div_gu), trunc_decim));

% ---- exclude the saddle point in rantzer's example 2 ----
% deg_safe_p8 = polynomialDegree(geo.poly_X_r);
% deg_safe_s8 = total_deg - deg_safe_p8;
% [sos_safety_prog, poly_safe_s8] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s8', Psi, deg_safe_s8,'true');
% poly_Xexcluded_2 = 0.16 - (x1-1.5)^2 - x2^2;

ineq_safe_2 = num_div_rhof_gu + poly_safe_s3 * geometry.poly_X_r - poly_safe_s2 * geometry.poly_X;% + poly_safe_s8 * poly_Xexcluded_2;
% ineq_safe_2 = div_pf_sym * ax_sol_sym^2 + poly_safe_s3 * domain_geo.poly_X_r;
% ineq_safe_2 = vpa(ineq_safe_2, 5)
% ineq_safe_2 = div_pf_sym * ax_sol_sym^2 + num_div_gu;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_2);

deg_safe_s7 =  deg_div_rhof_gu - deg_pXinit;
deg_diff = total_deg - deg_safe_s7
[sos_safety_prog, poly_safe_s7] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s7', Psi, deg_safe_s7,sparse);

ineq_safe_5 = num_div_rhof_gu - poly_safe_s7*geometry.poly_X_init;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_5);

% 3. rho <= 0 in X_u + (X_d - X)
ineq_safe_3 = - poly_safe_rho_sym - poly_safe_s6 * geometry.poly_X_u;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_3);

ineq_safe_4 = - poly_safe_rho_sym - poly_safe_s4*geometry.poly_X_d + poly_safe_s5*geometry.poly_X;
sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_4);

% 4. rho >=0 in X_r
% ineq_safe_5 = poly_safe_rho_sym - poly_safe_s7*domain_geo.poly_X_r;
% sos_safety_prog = sosineq(sos_safety_prog, ineq_safe_5);

% abs_c_safety = mpvar('abs_c_safety', [Q_safe_rho, 1]);
% abs_c_safety_sym = p2s(abs_c_safety);
% sos_safety_prog = sosdecvar(sos_safety_prog, abs_c_safety_sym);
% for  i_safe = 1:Q_safe_rho
%     i_safe
%     sos_safety_prog = sosineq(sos_safety_prog, abs_c_safety_sym(i_safe) - c_safe_rho_sym(i_safe));
%     sos_safety_prog = sosineq(sos_safety_prog, abs_c_safety_sym(i_safe) + c_safe_rho_sym(i_safe));
% end
% 
% %% add a L1-norm objective for the safety rho polynomial
% sparse_obj = sum(c_safe_rho_sym(1:Q_safe_rho));
% sos_safety_prog = sossetobj(sos_safety_prog, sparse_obj);

%% solve the safety verifycation feasibility problem
disp('----- solving the safety eventuality problem -----')
option_safety.solver = 'sedumi';
% option_safety.params.tol = 1e-15;

% option_safety.params.stepdif = 1;
% option_safety.params.vplot = 1;
% option_safety.params.alg = 2;
% option_safety.params.theta = 1;
% option_safety.params.beta = 0.9;

% solve and show info
% [sos_safety_prog, sos_safety_info] = solve_show_info(sos_safety_prog, option_safety);
[sos_safety_prog, sos_safety_info] = sossolve(sos_safety_prog, option_safety);
print_sol_info(sos_safety_info)

rho_safe_sol = sosgetsol(sos_safety_prog, poly_safe_rho_sym);
% rho_safe_sol = p2s(clean_polynomial(s2p(rho_safe_sol), trunc_decim))
% rho_safe_sol = vpa(rho_safe_sol, truncation_precision);

if rho_safe_sol == 0
    disp('Trivial solution of safety Rho: constant 0')
end

%% manually verify the solved polynomials
% disp('------- manually check the sos constraints -------')
% optimal control
% ab_sol_sym = ax_sol_sym * poly_b_sym;
% bc_sol_sym = cx_sol_sym * poly_b_sym;
% 
% sol_term(1) = vec_div(F.*ax_sol_sym, x);
% sol_term(2) = vec_div(G*cx_sol_sym, x);
% sol_term(3) = vec_div(F.*ab_sol_sym, x);
% sol_term(4) = vec_div(G*bc_sol_sym, x);

% constraint_sym_sol = (1+Alph)*poly_b_sym*sum(sol_term(1:2)) - Alph*sum(sol_term(3:4));

%  res_sys_constraints = findsos(constraint_sym_sol);
%  res_ax = findsos(ax_sol_sym);
 
% %  safety program
%  safe_s1_sol = sosgetsol(sos_safety_prog, poly_safe_s1);
%  safe_s2_sol = sosgetsol(sos_safety_prog, poly_safe_s2);
%  safe_s3_sol = sosgetsol(sos_safety_prog, poly_safe_s3);
%  safe_s4_sol = sosgetsol(sos_safety_prog, poly_safe_s4);
%  safe_s5_sol = sosgetsol(sos_safety_prog, poly_safe_s5);
%  safe_s6_sol = sosgetsol(sos_safety_prog, poly_safe_s6);
%  
%  ineq_safe_1_sol = sosgetsol(sos_safety_prog, ineq_safe_1);
%  ineq_safe_2_sol = sosgetsol(sos_safety_prog, ineq_safe_2);
%  ineq_safe_3_sol = sosgetsol(sos_safety_prog, ineq_safe_3);
% %  ineq_safe_4_sol = sosgetsol(sos_safety_prog, ineq_safe_4);
%  
%  res_s1 = findsos(safe_s1_sol);
%  res_s2 = findsos(safe_s2_sol);
%  res_s3 = findsos(safe_s3_sol);
%  res_s4 = findsos(safe_s4_sol);
%  res_s5 = findsos(safe_s5_sol);
%  res_s6 = findsos(safe_s6_sol);
%  
%  res_safe_ineq1 = findsos(ineq_safe_1_sol);
%  res_safe_ineq2 = findsos(ineq_safe_2_sol);
%  res_safe_ineq3 = findsos(ineq_safe_3_sol);
% %  res_safe_ineq4 = findsos(ineq_safe_4_sol);
%  
%  if isempty(res_sys_constraints) || isempty(res_s1) || isempty(res_s2) || isempty(res_s3) || isempty(res_s4) || ...
%          isempty(res_s5) || isempty(res_s6) || isempty(res_safe_ineq1) || isempty(res_safe_ineq2) || isempty(res_safe_ineq3)% || isempty(res_safe_ineq4) %  || isempty(res_ax)
%  disp('Infeasible solution!')
%  keyboard
%  end

%  if isempty(res_safe_ineq1) || isempty(res_safe_ineq2) || isempty(res_safe_ineq3)
%      disp('Infeasible solution!')
%      keyboard
%  else
%      disp('all the sos constraints are satisfied')
%  end

 
%% saving data
load(ocp_data_file);
disp('----- Saving solved data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_file_name = ['experiments/',date,'_lbd_',num2str(lambda),'_ocp_feasratio_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.mat'];
disp(data_file_name)
save(data_file_name)

