function output = plot_ocp_results(filename)
%% load data
disp('======== function plot_result =========')

% load data
filename
load(filename)

geometry = domain_definition();

% elements to be plotted
y1i = plotting_elements.meshgrid.y1i;
y2i = plotting_elements.meshgrid.y2i;
ucont = plotting_elements.v_field.ucont;
vcont = plotting_elements.v_field.vcont;
xinits = plotting_elements.xints;
ycont = plotting_elements.trajectory.ycont;
ax_over_bx = plotting_elements.surfaces.ax_over_bx;
div_ocp = plotting_elements.surfaces.div_ocp;
h_0 = plotting_elements.surfaces.h_0;
% figure title
ttl = figure_title(dynamics_option, lambda, ocp_info, []);

%% fig 1 :trajectories and velocity field
fig1 = figure();
plot_sets(geometry);
title(ttl, 'Interpreter', 'none')

% vector field plotting
quiver(y1i, y2i, ucont, vcont);

% initial points plotting
for i_xinits = xinits(1:end)
    i_xinits = i_xinits(1);
    i_xinits = i_xinits{1};
    for i_point_x = 1:length(i_xinits)
        plot(i_xinits(i_point_x, 1), i_xinits(i_point_x, 2), 'r*')
    end
end

% plot trajectories
for i_ycont = ycont(1:end)
    i_ycont = i_ycont(1);
    i_ycont = i_ycont{1};
    cellfun(@(x) plot(x(:,1),x(:,2),'k'), i_ycont);
    xlim([-4.5, 4.5])
    ylim([-4.5, 4.5])
    % zlim([-10, 10]);
end
%% fig1 zoomed trajectories
fig_zoomed = figure();
plot_sets(geometry);
quiver(y1i, y2i, ucont, vcont);

% red mark for the initial points
for i_xinits = xinits(1:end)
    i_xinits = i_xinits(1);
    i_xinits = i_xinits{1};
    for i_point_x = 1:length(i_xinits)
        plot(i_xinits(i_point_x, 1), i_xinits(i_point_x, 2), 'r*')
    end
end

% plot trajectories
for i_ycont = ycont(1:end)
    i_ycont = i_ycont(1);
    i_ycont = i_ycont{1};
    cellfun(@(x) plot(x(:,1),x(:,2),'k'), i_ycont);
end

% -------- config vdp ---------
%exp1
% xlim([-1.5, 2.3])
% ylim([-2.3, 1.2])
% exp2
% xlim([-2.0, 1.6])
% ylim([-1.2, 2.0])
% zlim([-10, 10]);
% -------- integrator ---------
% xlim([-4.0, 1.0])
% ylim([-1.0, 4.0])
% zlim([-10, 10]);
% vdp-L1
xlim([-2.2, 1.0])
ylim([-1.0, 3.2])

%% fig 2: div(rho (f+gu)) surface plot for optimal control
fig2 = figure('visible','off');
% fig2 = figure();
plot_sets(geometry);
% title(ttl, 'Interpreter', 'none');

% div_LHS1_simp = vpa(div_LHS1, 10);
% fhdl = matlabFunction(div_ocp);
% div_rho_near_origin = feval(fhdl, 1e-10, 1e-10)

% surface plotting
hold on
fsurf(div_ocp)

xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
view(-19,56)

%% fig 3: h_0 surface for optimal control
fig3 = figure('visible','off');
% fig3 = figure();
plot_sets(geometry);
title(ttl, 'Interpreter', 'none');

% surface plotting
hold on
fsurf(h_0);

xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
view(-19,56)

%% fig 4:a/b surface for optimal control
fig4 = figure('visible','off');
plot_sets(geometry);
title(ttl, 'Interpreter', 'none');

% surface plotting
hold on
fsurf(ax_over_bx);

xlim([-4.5, 4.5])
ylim([-4.5, 4.5])
view(-19,56)

%% save optimal control results
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_driven_type = option.data_driven_option.type;
regularizing_type = option.control_penalty_type;
file_prefix = ['res/',date, '_', dynamics_option,'_trj_',regularizing_type,'_',data_driven_type,'_','lamb_',num2str(lambda),'ocp_feas_',num2str(ocp_info.feasratio)];
name_ocp_trj = [file_prefix,'.png'];
name_ocp_trj_fig = [file_prefix,'.fig'];
name_ocp_trj_zoomed = [file_prefix,'_zoomed','.png'];
name_ocp_trj_fig_zoomed = [file_prefix,'_zoomed','.fig'];

saveas(fig1, name_ocp_trj);
saveas(fig1, name_ocp_trj_fig);

saveas(fig_zoomed, name_ocp_trj_zoomed);
saveas(fig_zoomed, name_ocp_trj_fig_zoomed);

end

