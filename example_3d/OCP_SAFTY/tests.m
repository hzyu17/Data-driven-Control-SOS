% test polynomial semi-algebraic set
syms x1 x2
[X_0, X, X_d, X_init, X_u, X_r] = domain_definition();

prog = sosprogram([x1, x2]);

total_deg = 12;
nx = 2;
nPsi = nchoosek(total_deg+nx, nx);
Psi = monomials(x, 0:total_deg);

deg_x = 10;
Qx = nchoosek(deg_x+nx, nx);

pvar x1 x2
x = [x1;x2];

c_x = mpvar('c_x', [Qx, 1]);
c_x = [c_x(1:Qx); zeros(nPsi-Qx, 1)];
c_x_sym = p2s(c_x);
poly_x = c_x'* Psi;
poly_x_sym = p2s(poly_x);
prog = sosdecvar(prog, c_x_sym(1:Qx));

% try with rational
% poly_b_sym = 

eps = 1e-9;
deg_s1 = deg_x - 2;

Qs1 = nchoosek(total_deg+nx, nx);
c_s1 = mpvar('c_s1', [Qs1, 1]);
c_s1 = [c_s1(1:Qs1); zeros(nPsi-Qs1,1)];
s1 = p2s(c_s1' * Psi);
c_s1_sym = p2s(c_s1);
prog = sosdecvar(prog, c_s1_sym(1:Qs1));

c_s2 = mpvar('c_s2', [Qs1, 1]);
c_s2 = [c_s2(1:Qs1); zeros(nPsi-Qs1,1)];
s2 = p2s(c_s2' * Psi);
c_s2_sym = p2s(c_s2);
% prog = sosdecvar(prog, c_s2_sym(1:Qs1));

c_s3 = mpvar('c_s3', [Qs1, 1]);
c_s3 = [c_s3(1:Qs1); zeros(nPsi-Qs1,1)];
s3 = p2s(c_s3' * Psi);
c_s3_sym = p2s(c_s3);
prog = sosdecvar(prog, c_s3_sym(1:Qs1));

c_s4 = mpvar('c_s4', [Qs1, 1]);
c_s4 = [c_s4(1:Qs1); zeros(nPsi-Qs1,1)];
s4 = p2s(c_s4' * Psi);
c_s4_sym = p2s(c_s4);
% prog = sosdecvar(prog, c_s4_sym(1:Qs1));

syms x1 x2
X_init = 4 - (x1+2)^2 - (x2-2)^2;
X_u = 4 - (x1-2)^2 - (x2+2)^2;
X_d = 49 - x1^2 - x2^2;

prog = sosineq(prog, s1);
% prog = sosineq(prog, s2);
prog = sosineq(prog, s3);
% prog = sosineq(prog, s4);

prog = sosineq(prog, poly_x_sym);
prog = sosineq(prog, poly_x_sym-0.9 - s1*X_init);
% prog = sosineq(prog, -poly_x_sym+1 - s2*X_init);
prog = sosineq(prog, -poly_x_sym+0.1 + s3*X_u);

% prog = sosineq(prog, poly_x_sym + s3*X_init);
% prog = sosineq(prog, poly_x_sym - eps - s4*X_u);
% prog = sosineq(prog, poly_x_sym);

obj = ones([Qx, 1])'*c_x_sym(1:Qx);
prog = sossetobj(prog, obj);

[prog, info] = sossolve(prog);

cx_sol_sym = sosgetsol(prog, poly_x_sym)

%%
figure
grid on
hold on
fsurf(cx_sol_sym)
% xlim([-4, 0])
% ylim([0, 4])
% zlim([0, 1])
xlabel('x1')
ylabel('x2')

figure
grid on
hold on
fsurf(cx_sol_sym)
% xlim([-4, 0])
% ylim([0, 4])
zlim([0, 1])
xlabel('x1')
ylabel('x2')

%% test on integration inside a circle
syms x1 x2
[X_0, X, X_d, X_init, X_u, X_r] = domain_definition();

%%
lambda = 50;
eps = 0.1;
rat_poly = 'Rational'
plot_result(lambda, eps, rat_poly)




