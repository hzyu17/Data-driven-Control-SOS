function [outputArg1,outputArg2] = plot_results(lambda, eps, poly_flg, geo)
%% load data
% load('sos_prog.mat')
disp('Loading solved data to plot')
if poly_flg == 1
    name = ['experiments/results_SOS_Poly_lbd_',num2str(lambda),'_eps_',num2str(eps),'.mat'];
else
    name = ['experiments/results_SOS_Rational_lbd_',num2str(lambda),'_eps_',num2str(eps),'.mat'];
end
load(name)

X_init_cen = geo.X_init_cen;
X_u_cen = geo.X_u_cen;
X_r_cen = geo.X_r_cen;
X_d_cen = geo.X_d_cen;
X_cen = geo.X_cen;

r_xinit = geo.r_xinit;
r_xr = geo.r_xr;
r_xu = geo.r_xu;
r_X = geo.r_X;
r_Xd = geo.r_Xd;

%% calculations
syms x1 x2
% ============
% Draw a vector field 
% ============
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
    xdot = feval(@(pt) control_dynamics_VDP(0, pt, x, cx_sol, ax_sol), yi(i1,: )');
    ucont(i1)= xdot(1);
    vcont(i1) = xdot(2);
end

% -------------------------
% plot result trajectories 
% -------------------------
pvar x1 x2
x = [x1; x2];
xinits = [ -0.4,-0.4; -0.2,-0.2; -0.0,0.3; 0.2,0.1; 0.3, 0.4; 0.6, 0.6; 1, 0; 0, 1; -1, 0]*sqrt(r_xinit) + X_init_cen;

% -----------------------------------------------------
% 2. solving ODE and stablizing using ux = cx / ax
% -----------------------------------------------------
tspan = [0 40];
for l1 = 1:length(xinits)
    xinit = xinits(l1, :);
    [tcont{l1}, ycont{l1}] = ode45(@(t,pt) control_dynamics_VDP(t, pt, x, cx_sol, ax_sol), tspan, xinit');
end

%% fig 1
% fig1 = figure('visible', 'off');
fig1 = figure()
ttl=['lambda', num2str(lambda), ' - eps', num2str(eps)];
title(ttl)
quiver(y1i, y2i, ucont, vcont);
xlabel('x_1');
ylabel('x_2');
axis tight equal;

hold on
grid on

% plot sets
syms x1 x2
x=[x1;x2];
X_init_cen_fun = @(x1, x2) r_xinit - (x1-X_init_cen(1)).^2 - (x2-X_init_cen(2)).^2;
fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
fc_init.LevelList = 0;
t_x0 = text(X_init_cen(1)+0.2, X_init_cen(2)+0.6,'X_0','Color','k','FontSize',15)
t_x0.FontWeight = 'bold';

X_u_cen_fun = @(x1, x2) r_xu - (x1-X_u_cen(1)).^2 - (x2-X_u_cen(2)).^2;
fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
fc_u.LevelList = 0;
t_xu = text(X_u_cen(1)-0.3, X_u_cen(2)-0.2, 'X_u','Color','r','FontSize',15)
t_xu.FontWeight = 'bold';

X_r_cen_fun = @(x1, x2) r_xr - (x1-X_r_cen(1)).^2 - (x2-X_r_cen(2)).^2;
fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
fc_r.LevelList = 0;
t_xr = text(X_r_cen(1)+0.4, X_r_cen(2)+0.5, 'X_r','Color','g','FontSize',15)
t_xr.FontWeight = 'bold';

X_d_cen_fun = @(x1, x2) r_Xd - (x1-X_d_cen(1)).^2 - (x2-X_d_cen(2)).^2;
fc_xd = fcontour(X_d_cen_fun, 'm','LineWidth',2);
fc_xd.LevelList = 0;
t_xd=text(X_d_cen(1)+3.0, X_d_cen(2)+3.0, 'X_d','Color','m','FontSize',16)
t_xd.FontWeight = 'bold';

X_cen_fun = @(x1, x2) r_X - (x1-X_cen(1)).^2 - (x2-X_cen(2)).^2;
fc_X = fcontour(X_cen_fun, 'b','LineWidth',2);
fc_X.LevelList = 0;
t_x=text(X_cen(1)+2.5, X_cen(2)+2.5, 'X','Color','b','FontSize',16)
t_x.FontWeight = 'bold';
legend('X')


for i_point_x = 1:length(xinits)
    plot(xinits(i_point_x, 1), xinits(i_point_x, 2), 'r*')
end

% plot trajectories
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
zlim([-1000, 1000]);
taglbd = ['Lambda = ', num2str(lambda)];
t_lbd=text(-1.0, 3.0, taglbd,'Color','k','FontSize',16);
t_lbd.FontWeight = 'bold';

if lambda >0
    hold on
    % Level set safety Rho
    Bx = matlabFunction(rho_safe_sol);
    fc = fcontour(rho_safe_sol,'--k',[-4 4],'LineWidth',2);
    fc.LevelList = 0;

    % level set of div(p(f+gu))
    div_all = vec_div(rho_safe_sol.*(F+G.*cx_sol_sym/ax_sol_sym), [x1, x2]);
    Bx1 = matlabFunction(div_all);
    fc_divall = fcontour(div_all,'--r',[-4 4],'LineWidth',2);
    fc_divall.LevelList = 0;
end

%% fig 2
fig2 = figure('visible', 'off');
ttl=['div(rho)--', 'lambda', num2str(lambda), ' - eps', num2str(eps)];
title(ttl);

if poly_flg == 1
    Rho = ax_sol_sym;
    Rho_bar = cx_sol_sym;
    LHS = Rho*F + Rho_bar*G;
    div_LHS1 = diff(LHS(1), x1) + diff(LHS(2), x2);
else % rational
    Rho = ax_sol_sym / ((bx_sol_sym)^Alph);
    Rho_bar = cx_sol_sym / ((bx_sol_sym)^Alph);
    LHS = Rho*F + Rho_bar*G;
    div_LHS1 = diff(LHS(1), x1) + diff(LHS(2), x2);
end

div_LHS1_simp = vpa(div_LHS1, 5);
fhdl = matlabFunction(div_LHS1_simp);
div_rho_near_origin = feval(fhdl, 1e-7, 1e-7)

hold on
fsurf(div_LHS1)

fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
fc_init.LevelList = 0;
text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14)

fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
fc_u.LevelList = 0;
text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14)

fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
fc_r.LevelList = 0;
text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14)

taglbd = ['Lambda = ', num2str(lambda)];
text(0, 4.0, taglbd,'Color','k','FontSize',16,'fontweight', 'bold');

xlabel('x1');
ylabel('x2');
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])

%% fig 3
fig3 = figure('visible','off');
ttl=['rho--', 'lambda', num2str(lambda), ' - eps', num2str(eps)];
title(ttl);
hold on
fsurf(Rho);

fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
fc_init.LevelList = 0;
text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14)

fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
fc_u.LevelList = 0;
text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14)

fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
fc_r.LevelList = 0;
text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14)

taglbd = ['Lambda = ', num2str(lambda)];
text(0, 4.0, taglbd,'Color','k','FontSize',16);

xlabel('x1');
ylabel('x2');
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])

%% fig4
% fig4 = figure('visible', 'off');
fig4 = figure();
fsurf(rho_safe_sol)
% zlim([-1, 1])
hold on
fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
fc_init.LevelList = 0;
text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14)

fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
fc_u.LevelList = 0;
text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14)

fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
fc_r.LevelList = 0;
text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14)

taglbd = ['Lambda = ', num2str(lambda)];
text(0, 4.0, taglbd,'Color','k','FontSize',16);

zlim([-1000, 1000]);

%% fig5
% div(rho*F + gu)
if lambda > 0
    fig5 = figure();
    fsurf(div_all)
    % zlim([-1, 1])
    hold on
    fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
    fc_init.LevelList = 0;
    text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14)

    fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
    fc_u.LevelList = 0;
    text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14)

    fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
    fc_r.LevelList = 0;
    text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14)

    taglbd = ['Lambda = ', num2str(lambda)];
    text(0, 4.0, taglbd,'Color','k','FontSize',16);

    zlim([-1000, 1000]);
end

%% save images
if poly_flg == 1
    name_trj = ['res/Poly_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.png'];
    name_trj_fig = ['res/Poly_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_div = ['res/Poly_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_rho = ['res/Poly_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_rho = ['res/Poly_safe_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
else
    name_trj = ['res/Rational_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.png'];
    name_trj_fig = ['res/Rational_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_div = ['res/Rational_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_rho = ['res/Rational_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_rho = ['res/Rational_safe_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_div = ['res/Rational_safe_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
end
saveas(fig1, name_trj);
saveas(fig1, name_trj_fig);
saveas(fig2, name_div);
saveas(fig3, name_rho);
saveas(fig4, name_safety_rho);
if lambda > 0
    saveas(fig5, name_safety_div);
end

finish = 'image saved'

end

