function ocp_plot_elements = plotting_handler(data_file,  experiment_option)
% calculate the elements to plot, e.g., the vector fields, the
% trajectories, and the surfaces, etc.

%%
if ~isempty(data_file)
    disp('Loading solved data to plot')
    data_file
    load(data_file)
else
    keyboard
end

% geometry = domain_definition();

if exist('var_poly') == 0
    pvar x1 x2
    var_poly.x = [x1; x2];
end

if exist('var_sym') == 0
    syms x1 x2
    var_sym.x = [x1; x2];
end

if exist('dynamics_sym') == 0
    dynamics_sym = dynamics_definition('sym', dynamics_option);
end
geometry = domain_definition();
%%
% ==========
% vector field 
% ==========
yl = -4.5;
yr = 4.5;
dy = 0.5;
[y1i, y2i] = meshgrid(yl:dy:yr, yl:dy:yr);
yi = [y1i(:), y2i(:)];
ucont = zeros(size(y1i));
vcont = zeros(size(y2i));

for i1 = 1:numel(y1i)
    xdot = feval(@(pt) control_dynamics_defn(dynamics_option, pt, var_poly.x, cx_sol, ax_sol), yi(i1,: )');
    ucont(i1) = xdot(1);
    vcont(i1) = xdot(2);
end

% ------ output ------
ocp_plot_elements.meshgrid.y1i = y1i;
ocp_plot_elements.meshgrid.y2i = y2i;
ocp_plot_elements.v_field.ucont = ucont;
ocp_plot_elements.v_field.vcont = vcont;

%%
% ---------------
% trajectories 
% ---------------
num_xinit = geometry.num_xinit;
ocp_plot_elements.xints = cell(1, num_xinit);
ocp_plot_elements.trajectory.tcont = cell(1, num_xinit);
ocp_plot_elements.trajectory.ycont = cell(1, num_xinit);
switch_controller = experiment_option.switch_controller;
for i_xinit = 1:num_xinit
    xinits = [ -0.4,-0.4; -0.2,-0.2; -0.0,0.3; 0.2,0.1; 0.3, 0.4; 0.6, 0.6; 1, 0; 0, 1; -1, 0]*geometry.r_xinit(i_xinit) + geometry.X_init_cen(i_xinit, :);
    ocp_plot_elements.xints{i_xinit} = xinits;

    % ----------------------------------------
    % 2. solving ODE using ux = cx / ax
    % ----------------------------------------
    tspan = [0 30];
    for l1 = 1:length(xinits)
        xinit = xinits(l1, :);
        [tcont{l1}, ycont{l1}] = ode15s(@(t,pt) control_dynamics_defn(dynamics_option, pt, var_poly.x, cx_sol, ax_sol), tspan, xinit');
    end
    %% local lqr controller
    if strcmp(switch_controller, 'local_lqr')
        for l1 = 1:length(xinits)
            start_point = 1;
            while sum(ycont{1, l1}(start_point, :).*ycont{1, l1}(start_point, :)) > geometry.r_xr^2 && start_point < size(ycont{1, l1}, 1)
                start_point = start_point + 1;
            end
            if start_point == size(ycont{1, l1}, 1) % no trajectory inside the X_r
                disp('no trajectory entering the reach set')
            end
            xinit = ycont{1, l1}(start_point, :);
            [~, ycont_lqr{l1}] = ode15s(@(t,pt) local_lqr_control(lc_lqr, pt, var_poly.x), tspan, xinit');
            ycont{l1} = [ycont{l1}(1:start_point, :); ycont_lqr{l1}];
        end
    end
    ocp_plot_elements.trajectory.tcont{i_xinit} = tcont;
    ocp_plot_elements.trajectory.ycont{i_xinit} = ycont;
end

%%
ax_over_bx = ax_sol_sym / ((poly_b_sym)^Alph);
Rho_bar = cx_sol_sym ./ ax_sol_sym;
LHS = ax_over_bx.*dynamics_sym.F + dynamics_sym.G*Rho_bar;
div_ocp = vec_div(LHS, var_sym.x);
h_0 = div_ocp / (poly_b_sym^Alph);

ocp_plot_elements.surfaces.ax_over_bx = ax_over_bx;
ocp_plot_elements.surfaces.div_ocp = div_ocp;
ocp_plot_elements.surfaces.h_0 = h_0;
% save('ocp_plot_elements.mat', 'ocp_plot_elements');

end

