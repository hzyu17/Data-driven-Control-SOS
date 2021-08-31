function data_file_name = create_safety_program()
%% sos program for safety verification
load('polynomials_def.mat')
geometry = domain_definition();

% define sos program
sos_safety_prog = sosprogram(var_sym.x);

% degrees
% VDP
% deg_safe_rho = 11;
deg_safe_rho = 11;

% deg_safe_rho = 3;
Q_safe_rho = nchoosek(deg_safe_rho+nx, nx);
c_safe_rho = mpvar('c_safe_rho', [Q_safe_rho, 1]);
c_safe_rho = [c_safe_rho(1:Q_safe_rho); zeros(nPsi-Q_safe_rho, 1)];
c_safe_rho_sym = p2s(c_safe_rho);
poly_safe_rho_sym = p2s(transpose(c_safe_rho) * Psi);

sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_rho_sym(1:Q_safe_rho));
sparse = 'false';

deg_pX = polynomialDegree(geometry.poly_X);
deg_safe_s1 = deg_safe_rho - deg_pX;
[sos_safety_prog, poly_safe_s1] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s1', Psi, deg_safe_s1,sparse);

% 2-dimensional case
deg_diff = total_deg - (deg_a + deg_c + deg_safe_rho - 1);

deg_safe_s2 = deg_safe_rho - deg_pX;
[sos_safety_prog, poly_safe_s2] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s2', Psi, deg_safe_s2,sparse);

deg_pXr = polynomialDegree(geometry.poly_X_r);
deg_safe_s3 = total_deg - deg_pXr;
[sos_safety_prog, poly_safe_s3] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s3', Psi, deg_safe_s3,sparse);

deg_pXd = polynomialDegree(geometry.poly_X_d);
deg_safe_s4 = deg_safe_rho - deg_pXd;
[sos_safety_prog, poly_safe_s4] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s4', Psi, deg_safe_s4,sparse);

deg_safe_p5 = polynomialDegree(geometry.poly_X);
deg_safe_s5 = deg_safe_rho
[sos_safety_prog, poly_safe_s5] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s5', Psi, deg_safe_s5,sparse);

deg_safe_s6 = deg_safe_rho;
[sos_safety_prog, poly_safe_s6] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s6', Psi, deg_safe_s6,sparse);

% deg_p7 = polynomialDegree(geo.poly_X_r);
% deg_safe_s7 = deg_safe_rho - deg_p7;
% [sos_safety_prog, poly_safe_s7] = declare_sos_var(sos_safety_prog, nx, 'c_safe_s7', Psi, deg_safe_s7,'true');


 %% saving data
disp('----- Saving solved optimal control data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_file_name = ['sos_program_data/',date,'_safety_verification_sos_prog.mat'];
disp(data_file_name)
save(data_file_name)
% save(data_file_name, 'sos_safety_prog', 'poly_safe_rho_sym', 'c_safe_rho_sym', 'c_safe_rho', 'Q_safe_rho', ...
%                                                'poly_safe_s1', 'poly_safe_s2', 'poly_safe_s3','poly_safe_s6')
end

