function [outputArg1,outputArg2] = plot_results(lambda, eps, poly_flg, geo, data_file_name)
%% load data
% load('sos_prog.mat')
disp('Loading solved data to plot')
data_file_name
load(data_file_name)

% CDC paper data
% exp1
% load('03_25_16_56_results_SOS_Rational_lbd_0_eps_0.mat')
% exp2
% load('experiments/03_25_17_29_results_SOS_Rational_lbd_100000_eps_0.mat')

X_init_cen = geo.X_init_cen;
X_u_cen = geo.X_u_cen;
X_r_cen = geo.X_r_cen;
X_d_cen = geo.X_d_cen;
X_cen = geo.X_cen;

r_xinit = geo.r_xinit;

poly_X_init =geo.poly_X_init;
poly_X_u = geo.poly_X_u;
poly_X_r = geo.poly_X_r;
poly_X_d = geo.poly_X_d;
poly_X = geo.poly_X;

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
sin_x1 = load('sin_x1.mat');
sin_x = sin_x1.sin_x1;
for i1 = 1:numel(y1i)
    xdot = feval(@(pt) control_dynamics_VDP(0, pt, x, cx_sol, ax_sol, sin_x), yi(i1,: )');
    ucont(i1)= xdot(1);
    vcont(i1) = xdot(2);
end

% -------------------------
% plot result trajectories 
% -------------------------
pvar x1 x2
x = [x1; x2];
xinits = [ -0.4,-0.4; -0.2,-0.2; -0.0,0.3; 0.2,0.1; 0.3, 0.4; 0.6, 0.6; 1, 0; 0, 1; -1, 0]*r_xinit + X_init_cen;

% -----------------------------------------------------
% 2. solving ODE and stablizing using ux = cx / ax
% -----------------------------------------------------
tspan = [0 40];
for l1 = 1:length(xinits)
    xinit = xinits(l1, :);
    [tcont{l1}, ycont{l1}] = ode45(@(t,pt) control_dynamics_VDP(t, pt, x, cx_sol, ax_sol, sin_x), tspan, xinit');
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

%% plot sets
syms x1 x2
x=[x1;x2];
X_init_cen_fun = matlabFunction(poly_X_init);
fc_init = fcontour(X_init_cen_fun, 'c','LineWidth',3);
fc_init.LevelList = 0;
t_x0 = text(X_init_cen(1)+0.2, X_init_cen(2)-0.6,'X_0','Color','c','FontSize',15);
t_x0.FontWeight = 'bold';

% X_u_cen_fun = @(x1, x2) r_xu^2 - (x1-X_u_cen(1)).^2 - (x2-X_u_cen(2)).^2;

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

% plot trajectories
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
zlim([-10, 10]);
% taglbd = ['Lambda = ', num2str(lambda)];
% t_lbd=text(-1.0, 3.0, taglbd,'Color','k','FontSize',16);
% t_lbd.FontWeight = 'bold';

% if lambda >0
hold on
% Level set safety Rho
Bx = matlabFunction(rho_safe_sol);
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = 0;

% level set of div(p(f+gu))
div_all = vec_div(rho_safe_sol.*(F+G*cx_sol_sym./ax_sol_sym), [x1, x2]);
Bx1 = matlabFunction(div_all);
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = 0;
% end

%% zoomed
fig_1_zoomed = figure();
quiver(y1i, y2i, ucont, vcont);
xlabel('x_1');
ylabel('x_2');
axis tight equal;

hold on
grid on

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
t_xr = text(X_r_cen(1)-0.5, X_r_cen(2)-0.6, 'X_r','Color','b','FontSize',15);
t_xr.FontWeight = 'bold';

fc_xd = fcontour(X_d_cen_fun, 'm','LineWidth',2);
fc_xd.LevelList = 0;
t_xd=text(X_d_cen(1)+3.0, X_d_cen(2)+4.0, 'X_d','Color','m','FontSize',15);
t_xd.FontWeight = 'bold';

fc_X = fcontour(X_cen_fun,'-b', 'LineWidth',2);
fc_X.LevelList = 0;
t_x=text(X_cen(1)+2.5, X_cen(2)-2.5, 'X','Color',[0.39,0.83,0.07],'FontSize',15);
t_x.FontWeight = 'bold';

for i_point_x = 1:length(xinits)
    plot(xinits(i_point_x, 1), xinits(i_point_x, 2), 'r*')
end

% plot trajectories
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
% -------- vdp ---------
%exp1
% xlim([-1.5, 2.3])
% ylim([-2.3, 1.2])
% exp2
xlim([-2.0, 1.6])
ylim([-1.2, 2.0])
zlim([-10, 10]);
% -------- integrator ---------
xlim([-1.0, 4.0])
ylim([-1.0, 3.5])
zlim([-10, 10]);
% taglbd = ['Lambda = ', num2str(lambda)];
% t_lbd=text(-1.0, 3.0, taglbd,'Color','k','FontSize',16);
% t_lbd.FontWeight = 'bold';

% if lambda >0
hold on
% Level set safety Rho
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = 0;

% level set of div(p(f+gu))
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = 0;

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
    Rho_bar = cx_sol_sym ./ ((bx_sol_sym)^Alph);
    LHS = Rho*F + G*Rho_bar;
    div_LHS1 = diff(LHS(1), x1) + diff(LHS(2), x2);
end

div_LHS1_simp = vpa(div_LHS1, 7);
fhdl = matlabFunction(div_LHS1_simp);
div_rho_near_origin = feval(fhdl, 1e-7, 1e-7)

hold on
fsurf(div_LHS1)

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

% fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
% fc_init.LevelList = 0;
% text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14);
% 
% fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
% fc_u.LevelList = 0;
% text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14);
% 
% fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
% fc_r.LevelList = 0;
% text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14);

% taglbd = ['Lambda = ', num2str(lambda)];
% text(0, 4.0, taglbd,'Color','k','FontSize',16,'fontweight', 'bold');

xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])

%% fig 3
fig3 = figure('visible','off');
ttl=['rho--', 'lambda', num2str(lambda), ' - eps', num2str(eps)];
title(ttl);
hold on
fsurf(Rho);

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

% fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
% fc_init.LevelList = 0;
% text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14);
% 
% fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
% fc_u.LevelList = 0;
% text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14);
% 
% fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
% fc_r.LevelList = 0;
% text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14);

taglbd = ['Lambda = ', num2str(lambda)];
text(0, 4.0, taglbd,'Color','k','FontSize',16);

xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])

%% fig4
% fig4 = figure('visible', 'off');
fig4 = figure();
fsurf(rho_safe_sol, 'k', 'FaceAlpha', 0.3)
% zlim([-1, 1])
hold on

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

% fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
% fc_init.LevelList = 0;
% text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14);
% 
% fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
% fc_u.LevelList = 0;
% text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14);
% 
% fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
% fc_r.LevelList = 0;
% text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14);

% taglbd = ['Lambda = ', num2str(lambda)];
% text(0, 4.0, taglbd,'Color','k','FontSize',16);

xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
% -------------- vdp ---------------
% experiment 1
xlim([-1.3, 2.3])
ylim([-2.3, 1.3])
zlim([-2, 2]);

% experiment 2
% xlim([-2.2, 1.6])
% ylim([-0.8, 2.0])
% zlim([-2, 2]);

% ----------- integrator ------------
xlim([-1.0, 4.0])
ylim([-1.0, 3.5])
zlim([-10, 10]);

%% fig5
% div(rho*F + gu)
% if lambda > 0
    fig5 = figure();
    fsurf(div_all, 'y', 'FaceAlpha', 0.3)
    % zlim([-1, 1])
    hold on
    
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

%     fc_init = fcontour(X_init_cen_fun, 'k','LineWidth',2);
%     fc_init.LevelList = 0;
%     text(X_init_cen(1)+0.1, X_init_cen(2)+0.1,'X_0','Color','k','FontSize',14);
% 
%     fc_u = fcontour(X_u_cen_fun, 'r','LineWidth',2);
%     fc_u.LevelList = 0;
%     text(X_u_cen(1)+0.1, X_u_cen(2)+0.1, 'X_u','Color','r','FontSize',14);
% 
%     fc_r = fcontour(X_r_cen_fun, 'g','LineWidth',2);
%     fc_r.LevelList = 0;
%     text(X_r_cen(1)+0.1, X_r_cen(2)+0.1, 'X_r','Color','g','FontSize',14);
% 
%     taglbd = ['Lambda = ', num2str(lambda)];
%     text(0, 4.0, taglbd,'Color','k','FontSize',16);
    
    xlabel('x_1', 'FontSize', 16);
    ylabel('x_2', 'FontSize', 16);
    % experiment 1
    xlim([-1.3, 2.3])
    ylim([-2.3, 1.3])
    zlim([-2, 2]);
    % experiment 2
%     xlim([-2.2, 1.6])
%     ylim([-0.8, 2.0])
%     zlim([-3, 2]);

% end

%% save images
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
if poly_flg == 1
    name_trj = ['res/',date, '_Poly_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.png'];
    name_trj_fig = ['res/',date, '_Poly_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_div = ['res/',date, '_Poly_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_rho = ['res/',date, '_Poly_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_rho = ['res/',date, '_Poly_safe_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
else
    name_trj = ['res/',date, '_Rational_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.png'];
    name_trj_fig = ['res/',date, '_Rational_trj_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_div = ['res/',date, '_Rational_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_rho = ['res/',date, '_Rational_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_rho = ['res/',date, '_Rational_safe_rho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
    name_safety_div = ['res/',date, '_Rational_safe_divrho_','lamb_',num2str(lambda),'_eps_',num2str(eps),'.fig'];
end
saveas(fig1, name_trj);
saveas(fig1, name_trj_fig);
saveas(fig2, name_div);
saveas(fig3, name_rho);
saveas(fig4, name_safety_rho);
% if lambda > 0
    saveas(fig5, name_safety_div);
% end

finish = 'image saved'

end

