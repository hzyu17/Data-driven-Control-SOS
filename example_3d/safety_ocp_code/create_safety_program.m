function [outputArg1,outputArg2] = create_safety_program()
%% sos program for safety verification
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

% define a new sos programsyms x1 x2
syms x1 x2
sos_safety_prog = sosprogram([x1, x2]);

deg_safe_rho = 12;
Q_safe_rho = nchoosek(deg_safe_rho+nx, nx);
c_safe_rho = mpvar('c_safe_rho', [Q_safe_rho, 1]);
c_safe_rho = [c_safe_rho(1:Q_safe_rho); zeros(nPsi-Q_safe_rho, 1)];
c_safe_rho_sym = p2s(c_safe_rho);
poly_safe_rho_sym = p2s(c_safe_rho' * Psi);

sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_rho_sym(1:Q_safe_rho));

deg_safe_p1 = polynomialDegree(X);
deg_safe_s1 = deg_safe_rho - deg_safe_p1;
Q_safe_s1 = nchoosek(deg_safe_s1+nx, nx);
c_safe_s1 = mpvar('c_safe_s1', [Q_safe_s1, 1]);
c_safe_s1 = [c_safe_s1(1:Q_safe_s1); zeros(nPsi-Q_safe_s1,1)];
poly_safe_s1 = p2s(c_safe_s1' * Psi);
c_safe_s1_sym = p2s(c_safe_s1);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s1_sym(1:Q_safe_s1));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s1);

deg_safe_p2 = polynomialDegree(X);
deg_diff =  total_deg - (deg_safe_rho + 1 + 2*deg_a - deg_safe_p2)
deg_safe_s2 = total_deg - deg_safe_p2;
Q_safe_s2 = nchoosek(deg_safe_s2+nx, nx);
c_safe_s2 = mpvar('c_safe_s2', [Q_safe_s2, 1]);
c_safe_s2 = [c_safe_s2(1:Q_safe_s2); zeros(nPsi-Q_safe_s2, 1)];
poly_safe_s2 = p2s(c_safe_s2' * Psi);
c_safe_s2_sym = p2s(c_safe_s2);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s2_sym(1:Q_safe_s2));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s2);

deg_safe_p3 = polynomialDegree(X_r);
deg_safe_s3 = total_deg - deg_safe_p3;
Q_safe_s3 = nchoosek(deg_safe_s3+nx, nx);
c_safe_s3 = mpvar('c_safe_s3', [Q_safe_s3, 1]);
c_safe_s3 = [c_safe_s3(1:Q_safe_s3); zeros(nPsi-Q_safe_s3, 1)];
poly_safe_s3 = p2s(c_safe_s3' * Psi);
c_safe_s3_sym = p2s(c_safe_s3);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s3_sym(1:Q_safe_s3));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s3);

deg_safe_p4 = polynomialDegree(X_d);
deg_safe_s4 = deg_safe_rho - deg_safe_p4;
Q_safe_s4 = nchoosek(deg_safe_s4+nx, nx);
c_safe_s4 = mpvar('c_safe_s4', [Q_safe_s4, 1]);
c_safe_s4 = [c_safe_s4(1:Q_safe_s4); zeros(nPsi-Q_safe_s4, 1)];
poly_safe_s4 = p2s(c_safe_s4' * Psi);
c_safe_s4_sym = p2s(c_safe_s4);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s4_sym(1:Q_safe_s4));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s4);

deg_safe_p5 = polynomialDegree(X);
deg_safe_s5 = deg_safe_rho - deg_safe_p5;
Q_safe_s5 = nchoosek(deg_safe_s5+nx, nx);
c_safe_s5 = mpvar('c_safe_s5', [Q_safe_s5, 1]);
c_safe_s5 = [c_safe_s5(1:Q_safe_s5); zeros(nPsi-Q_safe_s5, 1)];
poly_safe_s5 = p2s(c_safe_s5' * Psi);
c_safe_s5_sym = p2s(c_safe_s5);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s5_sym(1:Q_safe_s5));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s5);

deg_safe_p6 = polynomialDegree(X_u);
deg_safe_s6 = deg_safe_rho - deg_safe_p6;
Q_safe_s6 = nchoosek(deg_safe_s6+nx, nx);
c_safe_s6 = mpvar('c_safe_s6', [Q_safe_s6, 1]);
c_safe_s6 = [c_safe_s6(1:Q_safe_s6); zeros(nPsi-Q_safe_s6, 1)];
poly_safe_s6 = p2s(c_safe_s6' * Psi);
c_safe_s6_sym = p2s(c_safe_s6);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s6_sym(1:Q_safe_s6));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s6);

deg_safe_p7 = polynomialDegree(X_r);
deg_safe_s7 = deg_safe_rho - deg_safe_p7;
Q_safe_s7 = nchoosek(deg_safe_s7+nx, nx);
c_safe_s7 = mpvar('c_safe_s7', [Q_safe_s7, 1]);
c_safe_s7 = [c_safe_s7(1:Q_safe_s7); zeros(nPsi-Q_safe_s7, 1)];
poly_safe_s7 = p2s(c_safe_s7' * Psi);
c_safe_s7_sym = p2s(c_safe_s7);
sos_safety_prog = sosdecvar(sos_safety_prog, c_safe_s7_sym(1:Q_safe_s7));
sos_safety_prog = sosineq(sos_safety_prog, poly_safe_s7);
save('safety_prog')
end

