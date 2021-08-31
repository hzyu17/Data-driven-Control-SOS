function output = plot_from_file(data_file)
%% load data
% load('sos_prog.mat')
disp('\n======== function plot_result =========')
disp('Loading solved data to plot')
data_file
load(data_file)

geometry = domain_definition();

if exist('plotting_elements') == 0
    plotting_elements = plotting_handler(data_file);
end

if exist('var_poly') == 0
    pvar x1 x2
    var_poly.x = [x1; x2];
end

if exist('var_sym') == 0
    syms x1 x2
    var_sym.x = [x1; x2];
end

save(data_file)

y1i = plotting_elements.meshgrid.y1i;
y2i = plotting_elements.meshgrid.y2i;
ucont = plotting_elements.v_field.ucont;
vcont = plotting_elements.v_field.vcont;
xinits = plotting_elements.xints;
ycont = plotting_elements.trajectory.ycont;

%% parameters
% the levelset plotting value for safety density rho_safety and the
% divergence: div(rho_safety*(f + gu))
% level_rho_safety = -1e-2;
% level_div_rho_safe = 1e-3;
level_rho_safety = 0;
level_div_rho_safe = 0;

% title: the numerical condition of the sos programs
% ttl = ['lambda:', num2str(lambda),'-ocp feasratio:',num2str(ocp_info.feasratio),'-safety feasratio:',num2str(sos_safety_info.feasratio)];
ttl = figure_title(dynamics_option, lambda, ocp_info, sos_safety_info);
%%
% CDC paper data
% exp1
% load('03_25_16_56_results_SOS_Rational_lbd_0_eps_0.mat')
% exp2
% load('experiments/03_25_17_29_results_SOS_Rational_lbd_100000_eps_0.mat')

% load('experiments/06_18_15_58_results_SOS_Rational_lbd_50000.mat')

%% fig 1 :trajectories and velocity field
fig1 = figure();
title(ttl, 'Interpreter', 'none');
plot_sets(geometry);
quiver(y1i, y2i, ucont, vcont);

for i_point_x = 1:length(xinits)
    plot(xinits(i_point_x, 1), xinits(i_point_x, 2), 'r*')
end

% plot trajectories
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
hold on
% Level set safety Rho
% rantzer's paper
% rho_safe_sol = 7.6152 - 12.1597 * x1 - 3.85628 * x2 - 18.5818 * x1^2 - 3.92213 * x1 * x2 - 17.2032*x2^2 ...
%     + 51.0454*x1^3 + 15.8062*x1^2*x2 + 45.3684*x1*x2^2 - 32.2778*x1^4 - 9.27854*x1^3*x2 - ...
%     28.5075*x1^2*x2^2 + 1.03161*x1^2*x2^3 - 3.64267*x2^5 + 5.75452*x1^6 - 0.498478*x1^4*x2^2 - ...
%     1.60486*x1^2*x2^4 - 6.96206*x2^6 - 1.23579*x1^7 + 2.33577*x1^5*x2^2 - 1.86292*x1^4*x2^3 + 0.503928*x1^7*x2 + 2.19945*x1^3*x2^5;

Bx = matlabFunction(rho_safe_sol);
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = level_rho_safety;

% level set of div(p(f+gu))
div_all = vec_div(rho_safe_sol.*(F+G*cx_sol_sym./ax_sol_sym), var_sym.x);
Bx1 = matlabFunction(div_all);
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = level_div_rho_safe;

% div_pf_sym = vec_div(rho_safe_sol.*F, [x1, x2]);
% num_div_gu = vec_div(G*cx_sol_sym.*rho_safe_sol, [x1;x2])*ax_sol_sym - transpose((G*cx_sol_sym.*rho_safe_sol)) * vec_grad(ax_sol_sym, [x1, x2]);
% div_numerator = div_pf_sym .* ax_sol_sym^2 + num_div_gu;
% fc_div_numerator = fcontour(div_numerator,'-.y',[-4 4],'LineWidth',2);
% fc_div_numerator.LevelList = 0;

xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
zlim([-10, 10]);

%% fig1 zoomed trajectories
% fig_1_zoomed = plot_sets(geometry, 'visible');
fig_1_zoomed = figure();
title(ttl, 'Interpreter', 'none');
plot_sets(geometry);
quiver(y1i, y2i, ucont, vcont);
for i_point_x = 1:length(xinits)
    plot(xinits(i_point_x, 1), xinits(i_point_x, 2), 'r*')
end

% plot trajectories
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);

%% safety rho and div_rho 0-contour
hold on
% Level set safety Rho
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = level_rho_safety;

% level set of div(p(f+gu))
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = level_div_rho_safe;

% -------- vdp ---------
%exp1
% xlim([-1.5, 2.3])
% ylim([-2.3, 1.2])
% exp2
% xlim([-2.0, 1.6])
% ylim([-1.2, 2.0])
% zlim([-10, 10]);
% -------- integrator ---------
xlim([-4.0, 1.0])
ylim([-1.0, 4.0])
% zlim([-10, 10]);

%% fig 2: div(rho (f+gu)) surface plot for optimal control
% fig2 = figure('visible', 'off');
% title(ttl);
% fig2 = plot_sets(geometry);
% 
% % if poly_flg == 1
% %     Rho = ax_sol_sym;
% %     Rho_bar = cx_sol_sym;
% %     LHS = Rho*F + Rho_bar*G;
% %     div_LHS1 = diff(LHS(1), x1) + diff(LHS(2), x2);
% % else % rational
% ax_over_bx = ax_sol_sym / ((bx_sol_sym)^Alph);
% Rho_bar = cx_sol_sym ./ ((bx_sol_sym)^Alph);
% LHS = ax_over_bx*F + G*Rho_bar;
% div_LHS1 = vec_div(LHS, [x1; x2]);
% % end
% 
% % div_LHS1_simp = vpa(div_LHS1, 10);
% fhdl = matlabFunction(div_LHS1);
% div_rho_near_origin = feval(fhdl, 1e-10, 1e-10)
% 
% hold on
% fsurf(div_LHS1)
% 
% % xlabel('x1', 'FontSize', 14);
% % ylabel('x2', 'FontSize', 14);
% xlim([-4.5, 4.5])
% ylim([-4.5, 4.5])

%% fig 3: rho(f+gu) surface for optimal control
% fig3 = figure('visible', 'off');
% title(ttl);
% fig3 = plot_sets(geometry);
% hold on
% fsurf(ax_over_bx);
% 
% taglbd = ['Lambda = ', num2str(lambda)];
% text(0, 4.0, taglbd,'Color','k','FontSize',16);
% 
% % xlabel('x1', 'FontSize', 14);
% % ylabel('x2', 'FontSize', 14);
% xlim([-4.5, 4.5])
% ylim([-4.5, 4.5])

%% fig4: rho_safety surface plot
fig4 = figure();
title(ttl, 'Interpreter', 'none');
plot_sets(geometry);

% vector fields
hold on
quiver(y1i, y2i, ucont, vcont);

% surface plot
fsurf(rho_safe_sol, 'k', 'FaceAlpha', 0.3)

% surface plot
fsurf(div_all, 'y', 'FaceAlpha', 0.3)

% -------------- vdp ---------------
% experiment 1
% xlim([-1.3, 2.3])
% ylim([-2.3, 1.3])
% zlim([-2, 2]);
% experiment 2
% xlim([-2.2, 1.6])
% ylim([-0.8, 2.0])
% zlim([-2, 2]);
% ----------- integrator ------------
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
zlim([-10, 10]);
view(-19,56)

%% fig4: rho_safety surface plot: no z-limit
fig8 = figure();
title(ttl, 'Interpreter', 'none');
plot_sets(geometry);

% vector fields
hold on
quiver(y1i, y2i, ucont, vcont);

% surface plot
fsurf(rho_safe_sol, 'k', 'FaceAlpha', 0.3)

% surface plot
fsurf(div_all, 'y', 'FaceAlpha', 0.3)

% -------------- vdp ---------------
% experiment 1
% xlim([-1.3, 2.3])
% ylim([-2.3, 1.3])
% zlim([-2, 2]);
% experiment 2
% xlim([-2.2, 1.6])
% ylim([-0.8, 2.0])
% zlim([-2, 2]);
% ----------- integrator ------------
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
% zlim([-10, 10]);
view(-19,56)

%% a combined subplots for trajectory and safety verification
fig6 = figure('Position', [100, 100, 2048, 1200]);
s1 = subplot(2,2,[1 3]);
title(ttl, 'Interpreter', 'none');
xlim([-4.5, 4.5])
ylim([-4.5, 4.5])

plot_sets(geometry);
grid minor

hold on
quiver(y1i, y2i, ucont, vcont);
cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);

% Level set safety Rho
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = level_rho_safety;

% level set of div(p(f+gu))
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = level_div_rho_safe;

s2 = subplot(2,2,2);
view(-19,56)
plot_sets(geometry);
grid minor
hold on
quiver(y1i, y2i, ucont, vcont);
zlim([-10,10]);
fsurf(rho_safe_sol, 'k', 'FaceAlpha', 0.3)

s3 = subplot(2,2,4);
view(-19,56)
plot_sets(geometry);
grid minor
hold on
quiver(y1i, y2i, ucont, vcont);
zlim([-10,10]);
fsurf(div_all, 'y', 'FaceAlpha', 0.3)

%% trajectory, safety verifications in one figure
fig7 = figure('Position', [100, 100, 2048, 1200]);
title(ttl, 'Interpreter', 'none');
plot_sets(geometry);
grid minor
hold on
quiver(y1i, y2i, ucont, vcont);

cellfun(@(x) plot(x(:,1),x(:,2),'k'), ycont);
% Level set safety Rho
fc = fcontour(rho_safe_sol,'-.k',[-4 4],'LineWidth',2);
fc.LevelList = level_rho_safety;

% level set of div(p(f+gu))
fc_divall = fcontour(div_all,'-.y',[-4 4],'LineWidth',2);
fc_divall.LevelList = level_div_rho_safe;

view(-19,56)
zlim([-10,10]);
fsurf(rho_safe_sol, 'k', 'FaceAlpha', 0.1)
fsurf(div_all, 'y', 'FaceAlpha', 0.1)

%% save images
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
% if poly_flg == 1
%     name_trj = ['res/',date, '_Poly_trj_','lamb_',num2str(lambda),'.png'];
%     name_trj_fig = ['res/',date, '_Poly_trj_','lamb_',num2str(lambda),'.fig'];
%     name_div = ['res/',date, '_Poly_divrho_','lamb_',num2str(lambda),'.fig'];
%     name_rho = ['res/',date, '_Poly_rho_','lamb_',num2str(lambda),'.fig'];
%     name_safety_rho = ['res/',date, '_Poly_safe_rho_','lamb_',num2str(lambda),'.fig'];
% else
name_trj = ['res/',date, '_', dynamics_option,'_trj_','lamb_',num2str(lambda),'ocp_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.png'];
name_trj_fig = ['res/',date, '_', dynamics_option, '_trj_','lamb_',num2str(lambda),'ocp_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_div = ['res/',date, '_', dynamics_option, '_divrho_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_rho = ['res/',date, '_', dynamics_option, '_rho_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_safety_rho = ['res/',date, '_', dynamics_option, '_safe_rho_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_safety_rho_png = ['res/',date, '_', dynamics_option, '_safe_rho_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.png'];
name_safety_div = ['res/',date, '_', dynamics_option, '_safe_divrho_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_subplots_fig = ['res/',date, '_', dynamics_option, '_figures_combined_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_subplots_png = ['res/',date, '_', dynamics_option, '_figures_combined_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.jpg'];
name_stacked_fig = ['res/',date, '_', dynamics_option, '_figures_stacked_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.fig'];
name_stacked_png = ['res/',date, '_', dynamics_option, '_figures_stacked_','lamb_',num2str(lambda),'_feas_',num2str(ocp_info.feasratio),'_safety_feasratio_',num2str(sos_safety_info.feasratio),'.png'];
% end

saveas(fig1, name_trj);
saveas(fig1, name_trj_fig);
% saveas(fig2, name_div);
% saveas(fig3, name_rho);
saveas(fig4, name_safety_rho);
saveas(fig4, name_safety_rho_png);
% % saveas(fig5, name_safety_div);
saveas(fig6, name_subplots_fig);
saveas(fig6, name_subplots_png);
saveas(fig7, name_stacked_fig);
saveas(fig7, name_stacked_png);

disp('image saved')

end

