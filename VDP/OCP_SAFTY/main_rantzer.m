%% Test of Rantzer's method of generating a safe control for the trivial integrator dynamics
clear all
close all
clc

%% problem geometry definition
domain_geo = domain_definition();
poly_X_init = domain_geo.poly_X_init;
poly_X_u = domain_geo.poly_X_u;
poly_X_r = domain_geo.poly_X_r;
poly_X_d = domain_geo.poly_X_d;
poly_X = domain_geo.poly_X;
num_xu = domain_geo.num_xu;

X_init_cen = domain_geo.X_init_cen;
X_u_cen = domain_geo.X_u_cen;
X_r_cen = domain_geo.X_r_cen;
X_d_cen = domain_geo.X_d_cen;
X_cen = domain_geo.X_cen;
r_xinit = domain_geo.r_xinit;
r_xr = domain_geo.r_xr;

%% define polynomials
% VDP
total_deg = 14;
pvar x1 x2
x = [x1; x2];
nx = length(x);

% monomial
nPsi = nchoosek(total_deg+nx, total_deg);
Psi = monomials([x1 x2], [0:total_deg]);
Psi_sym = p2s(Psi);

% polynomial degrees
deg_a = 8;

% dimension of the control input c(x)
dim_m = 2;
deg_c = [10, 10];

% vdp
% dim_m = 1;
% deg_c = 6;

if length(deg_c) ~= dim_m
    disp('wrong dimension of inputs')
    keyboard
end

% ax 
Qa = nchoosek(deg_a+nx, nx);
c_a = mpvar('c_a', [Qa, 1]);
c_a = [c_a; zeros(nPsi - Qa, 1)];
c_a_sym = p2s(c_a);
% make sure that the constant in ax is 1
c_a_1 = sym(1);
c_a_sym = [c_a_1; c_a_sym(2:end)];

poly_a = transpose(c_a)* Psi;
poly_a_sym = transpose(c_a_sym) * Psi_sym;

% cx
Qc = [];
c_c = [];
c_c_sym = [];
poly_c = [];
poly_c_sym = [];

for i_c = 1:dim_m
    Qc_i = nchoosek(deg_c(i_c)+nx, nx);
    Qc = [Qc, Qc_i];
    var_name = ['c_c', num2str(i_c)];
    c_c_i = mpvar(var_name, [Qc_i, 1]);
    c_c_i = [c_c_i; zeros(nPsi - Qc_i, 1)];
    c_c = [c_c, {c_c_i}];
    c_c_i_sym = p2s(c_c_i);
    c_c_sym = [c_c_sym; {c_c_i_sym}];
    poly_c_i = transpose(c_c_i) * Psi;
    poly_c = [poly_c; poly_c_i];
    poly_c_sym_i = transpose(c_c_i_sym) * Psi_sym;
    poly_c_sym = [poly_c_sym; poly_c_sym_i];
end

%% sos program variable definition
geo = domain_definition();

% dynamics
[F, G] = dynamics_vdp();

if strcmp(class(F), 'polynomial')
    F = p2s(F);
    F = vpa(F, 8);
end

%% sos program decision variables declaration
syms x1 x2
x = [x1; x2];
sos_prog = sosprogram([x1, x2]);

for c_i = 1:dim_m
    c_c_sym_i = c_c_sym(c_i);
    c_c_sym_i = c_c_sym_i{1};
    Qci = Qc(c_i);
    sos_prog = sosdecvar(sos_prog, [c_c_sym_i(1:Qci)]);
end

c_a_sym = p2s(c_a);
sos_prog = sosdecvar(sos_prog, [c_a_sym(2:Qa)]);

deg_p1 = polynomialDegree(poly_X);
deg_safe_s1 = deg_a - deg_p1;
[sos_prog, poly_s1_sym] = declare_sos_var(sos_prog, nx, 'c_s1', Psi, deg_safe_s1);

poly_s2_syms = [];
for i_xu = 1:num_xu
    deg_p2_ixu = polynomialDegree(poly_X_u(i_xu));
    deg_s2_ixu = deg_a - deg_p2_ixu;
    var_name = ['c_s2_', num2str(i_xu)];
    [sos_prog, poly_s2_sym] = declare_sos_var(sos_prog, nx, var_name, Psi, deg_s2_ixu);
    poly_s2_syms = [poly_s2_syms, poly_s2_sym];
end

deg_safe_p3 = polynomialDegree(poly_X_d);
deg_safe_s3 = deg_a - deg_safe_p3;
[sos_prog, poly_s3_sym] = declare_sos_var(sos_prog, nx, 'c_s3', Psi, deg_safe_s3);

deg_safe_p4 = polynomialDegree(poly_X);
deg_safe_s4 = deg_a - deg_safe_p4;
[sos_prog, poly_s4_sym] = declare_sos_var(sos_prog, nx, 'c_s4', Psi, deg_safe_s4);

%% sos inequality construction and solving 
truncation_precision = 8;
trunc_decim = 1e-9;
eps = 1e3;
eps_div = 100;

[F, G] = dynamics_vdp();

ineq_bound_rho = eps - poly_a_sym;
sos_prog = sosineq(sos_prog, ineq_bound_rho);

% (1) div(rho(f+ug)) > 0 in X-X_r
term(1) = vec_div(F.*poly_a_sym, x);
term(2) = vec_div(G*poly_c_sym, x);

constraint_sym = sum(term(1:2));

ineq_bound_divrho = constraint_sym - eps_div;
sos_prog = sosineq(sos_prog, ineq_bound_divrho);

deg_safe_p5 = polynomialDegree(poly_X);
deg_constraint = polynomialDegree(constraint_sym)
deg_safe_s5 = deg_constraint - deg_safe_p5;
[sos_prog, poly_s5_sym] = declare_sos_var(sos_prog, nx, 'c_s5', Psi, deg_safe_s5);

deg_safe_p6 = polynomialDegree(poly_X_r);
deg_safe_s6 = polynomialDegree(constraint_sym) - deg_safe_p6;
[sos_prog, poly_s6_sym] = declare_sos_var(sos_prog, nx, 'c_s6', Psi, deg_safe_s6);

ineq_xr = constraint_sym - poly_s5_sym * poly_X + poly_s6_sym * poly_X_r;
% ineq_xr = p2s(clean_polynomial(s2p(ineq_xr), trunc_decim));
% ineq_xr = vpa(ineq_xr, truncation_precision);

sos_prog = sosineq(sos_prog, ineq_xr);
% sos_prog = sosineq(sos_prog, poly_a_sym);

% (2) rho > 0 in X0
ineq1 = poly_a_sym - poly_s1_sym * poly_X_init;
% ineq1 = p2s(clean_polynomial(s2p(ineq1), trunc_decim));
% ineq1 = vpa(ineq1, truncation_precision);
sos_prog = sosineq(sos_prog, ineq1);

% (3) rho <= 0 in Xu
for i_xu = 1:num_xu
    ineq2 = -poly_a_sym - poly_s2_syms(i_xu) * poly_X_u(i_xu);
    sos_prog = sosineq(sos_prog, ineq2);
end

% (4) rho <= 0 in Xd - X
ineq3 = -poly_a_sym - poly_s3_sym * poly_X_d + poly_s4_sym * poly_X;
sos_prog = sosineq(sos_prog, ineq3);

% solve the sos program
option.solver = 'sedumi';
option.params.tol = 1e-15;
sos_prog = solve_show_info(sos_prog, option);

%% obtain the solution
% get solution: cx, ax.
cx_sol = [];
cx_sol_sym = sosgetsol(sos_prog, poly_c_sym);
for i_c = 1:dim_m
    cx_sol_sym(i_c) = clean_polynomial(cx_sol_sym(i_c), trunc_decim);
%     cx_sol_sym(i_c) = vpa(cx_sol_sym(i_c), truncation_precision)
    if cx_sol_sym(i_c) == 0
        keyboard
    end
    cx_sol = [cx_sol; s2p(cx_sol_sym(i_c))]
end
ax_sol_sym = sosgetsol(sos_prog, poly_a_sym);
ax_sol_sym = clean_polynomial(s2p(ax_sol_sym), trunc_decim);
% ax_sol_sym = vpa(ax_sol_sym, truncation_precision)
ax_sol = s2p(ax_sol_sym)

%% plot the controlled system result
disp('============= plotting =============')

% calculations
syms x1 x2
% ===============
% Draw a vector field 
% ===============
yl = -4.5;
yr = 4.5;
dy = 0.5;
[y1i, y2i] = meshgrid(yl:dy:yr, yl:dy:yr);
yi = [y1i(:), y2i(:)];
ucont = zeros(size(y1i));
vcont = zeros(size(y2i));

pvar x1 x2
x = [x1; x2];
for i1 = 1:numel(y1i)
    xdot = feval(@(pt) control_dynamics_defn(pt, x, cx_sol, ax_sol), yi(i1,: )');
    ucont(i1)= xdot(1);
    vcont(i1) = xdot(2);
end

% ---------------------------------
%  trajectories starting points
% ---------------------------------
xinits = [ -0.4,-0.4; -0.2,-0.2; -0.0,0.3; 0.2,0.1; 0.3, 0.4; 0.6, 0.6; 1, 0; 0, 1; -1, 0]*r_xr + X_r_cen;
% xinits = [ -0.4,-0.4]*r_xinit + X_init_cen;

% --------------------------------------------------------
% 1. solving ODE and stablizing using ux = cx / ax
% --------------------------------------------------------
pvar x1 x2
x = [x1; x2];
disp('--------- solving odes ---------')
tspan = [0 40];
% cx_sol = x1 + x2; % test the ode45 function
% ax_sol  = 1; % test the ode45 function
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
for l1 = 1:length(xinits)
    xinit = xinits(l1, :);
    [tcont{l1}, ycont{l1}] = ode45(@(t,pt) control_dynamics_defn(pt, x, cx_sol, ax_sol), tspan, xinit')
end

%% fig 1
fig1 = figure();
ttl=['trajectory plot'];
title(ttl)
quiver(y1i, y2i, ucont, vcont);
xlabel('x_1');
ylabel('x_2');
axis tight equal;

hold on
grid on

%% plot sets
disp('--------- plotting sets ---------')
syms x1 x2
x=[x1;x2];
X_init_cen_fun = matlabFunction(poly_X_init);
fc_init = fcontour(X_init_cen_fun, 'c','LineWidth',3);
fc_init.LevelList = 0;
t_x0 = text(X_init_cen(1)+0.2, X_init_cen(2)-0.6,'X_0','Color','c','FontSize',15);
t_x0.FontWeight = 'bold';

X_u_cen_fun = matlabFunction(poly_X_u);
fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',3);
fc_u.LevelList = 0;
t_xu = text(X_u_cen(1)-0.3, X_u_cen(2)-0.2, 'X_u','Color','r','FontSize',15);
t_xu.FontWeight = 'bold';

X_r_cen_fun = matlabFunction(poly_X_r);
fc_r = fcontour(X_r_cen_fun, 'b','LineWidth',3);
fc_r.LevelList = 0;
t_xr = text(X_r_cen(1)-0.6, X_r_cen(2)-0.7, 'X_r','Color','b','FontSize',15);
t_xr.FontWeight = 'bold';

X_d_cen_fun = matlabFunction(poly_X_d);
fc_xd = fcontour(X_d_cen_fun, 'm','LineWidth',2);
fc_xd.LevelList = 0;
t_xd=text(X_d_cen(1)+3.0, X_d_cen(2)+4.0, 'X_d','Color','m','FontSize',15);
t_xd.FontWeight = 'bold';

X_cen_fun = matlabFunction(poly_X);
fc_X = fcontour(X_cen_fun,'-g', 'LineWidth',2);
fc_X.LevelList = 0;
t_x=text(X_cen(1)+2.5, X_cen(2)-2.5, 'X','Color',[0.39,0.83,0.07],'FontSize',15);
t_x.FontWeight = 'bold';

for i_point_x = 1:length(xinits)
    plot(xinits(i_point_x, 1), xinits(i_point_x, 2), 'r*')
end

%% plot trajectories
disp('---------- plotting trajectories ------------')
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
zlim([-10, 10]);

%% figure 2 
% rho
fig2 = figure();
title('rho')
fsurf(ax_sol_sym, 'b', 'FaceAlpha', 0.3)
% zlim([-1, 1])
hold on

fc_X = fcontour(X_cen_fun,'-g', 'LineWidth',2);
fc_X.LevelList = 0;
t_x=text(X_cen(1)+2.5, X_cen(2)-2.5, 'X','Color',[0.39,0.83,0.07],'FontSize',15);
t_x.FontWeight = 'bold';

fc_init = fcontour(X_init_cen_fun, 'c','LineWidth',3);
fc_init.LevelList = 0;
t_x0 = text(X_init_cen(1)+0.2, X_init_cen(2)+0.6,'X_0','Color','c','FontSize',15);
t_x0.FontWeight = 'bold';

fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',3);
fc_u.LevelList = 0;
t_xu = text(X_u_cen(1)-0.3, X_u_cen(2)-0.2, 'X_u','Color','r','FontSize',15);
t_xu.FontWeight = 'bold';

fc_r = fcontour(X_r_cen_fun, 'b','LineWidth',3);
fc_r.LevelList = 0;
t_xr = text(X_r_cen(1)+0.4, X_r_cen(2)+0.5, 'X_r','Color','b','FontSize',15);
t_xr.FontWeight = 'bold';

xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);

zlim([-1000, 1500]);

%% fig3
% div(rho*(F + gu))
term(1) = vec_div(F.*ax_sol_sym, x);
term(2) = vec_div(G*cx_sol_sym, x);

div_sol = sum(term(1:2));

fig3 = figure();
title('div(rho(f+gu))')
fsurf(div_sol, 'y', 'FaceAlpha', 0.3)
% zlim([-1, 1])
hold on

fc_X = fcontour(X_cen_fun,'-g', 'LineWidth',2);
fc_X.LevelList = 0;
t_x=text(X_cen(1)+2.5, X_cen(2)-2.5, 'X','Color',[0.39,0.83,0.07],'FontSize',15);
t_x.FontWeight = 'bold';

fc_init = fcontour(X_init_cen_fun, 'c','LineWidth',3);
fc_init.LevelList = 0;
t_x0 = text(X_init_cen(1)+0.2, X_init_cen(2)+0.6,'X_0','Color','c','FontSize',15);
t_x0.FontWeight = 'bold';

fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',3);
fc_u.LevelList = 0;
t_xu = text(X_u_cen(1)-0.3, X_u_cen(2)-0.2, 'X_u','Color','r','FontSize',15);
t_xu.FontWeight = 'bold';

fc_r = fcontour(X_r_cen_fun, 'b','LineWidth',3);
fc_r.LevelList = 0;
t_xr = text(X_r_cen(1)+0.4, X_r_cen(2)+0.5, 'X_r','Color','b','FontSize',15);
t_xr.FontWeight = 'bold';

xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);
% ---------------- vdp -------------------
% experiment 1
% xlim([-1.3, 2.3])
% ylim([-2.3, 1.3])
% zlim([-2, 2]);
% experiment 2
%     xlim([-2.2, 1.6])
%     ylim([-0.8, 2.0])
%     zlim([-3, 2]);
% --------------- integrator ---------------
% xlim([-4.0, 1.0])
% ylim([-1.0, 4.0])
zlim([-10, 1000]);

%% save images
disp('------------- saving images --------------')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
name_fig = ['res/rantzer/',date,'_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.fig'];
name_jpg = ['res/rantzer/',date,'_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.jpg'];
saveas(fig1, name_fig);
saveas(fig1, name_jpg);

name_fig = ['res/rantzer/',date,'_rho_3d_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.fig'];
name_jpg = ['res/rantzer/',date,'_rho_3d_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.jpg'];
saveas(fig2, name_fig);
saveas(fig2, name_jpg);

name_fig = ['res/rantzer/',date,'_divrho_3d_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.fig'];
name_jpg = ['res/rantzer/',date,'_divrho_3d_dega_',num2str(deg_a),'_degc_',num2str(deg_c(1)),'.jpg'];
saveas(fig3, name_fig);
saveas(fig3, name_jpg);
