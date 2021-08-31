function handle = plot_sets(geometry)
% Plot sets and experiment settings in the figure as a preprocessing.

% if strcmp(visible, 'visible')
%     fig = figure();
% elseif strcmp(visible, 'invisible')
%     fig = figure('visible', 'off');
% else
%     keyboard
% end
handle = gca;
num_xu = geometry.num_xu;
num_xinit = geometry.num_xinit;
% plot sets
xlabel('x1', 'FontSize', 16);
ylabel('x2', 'FontSize', 16);
% axis tight equal;

hold on
grid on

for i_xinit = 1:num_xinit
    X_init_cen_fun_i = matlabFunction(geometry.poly_X_init(i_xinit));
    fc_init_i = fcontour(X_init_cen_fun_i, 'c','LineWidth',3);
    fc_init_i.LevelList = 0;
    t_x0 = text(geometry.X_init_cen(i_xinit,1)+0.2, geometry.X_init_cen(i_xinit, 2)-0.6,'X_0','Color','c','FontSize',15);
    t_x0.FontWeight = 'bold';
end

for i_xu = 1:num_xu
    X_u_cen_fun_i = matlabFunction(geometry.poly_X_u(i_xu));
    fc_u_i = fcontour(X_u_cen_fun_i, 'r','LineWidth',3);
    fc_u_i.LevelList = 0;
    t_xu = text(geometry.X_u_cen(i_xu, 1)-0.3, geometry.X_u_cen(i_xu, 2)-0.2, 'X_u','Color','r','FontSize',15);
    t_xu.FontWeight = 'bold';
end

X_r_cen_fun = matlabFunction(geometry.poly_X_r);
fc_r = fcontour(X_r_cen_fun, 'b','LineWidth',3);
fc_r.LevelList = 0;
t_xr = text(geometry.X_r_cen(1)+0.4, geometry.X_r_cen(2)-0.4, 'X_r','Color','b','FontSize',15);
t_xr.FontWeight = 'bold';

X_d_cen_fun = matlabFunction(geometry.poly_X_d);
fc_xd = fcontour(X_d_cen_fun, 'm','LineWidth',2);
fc_xd.LevelList = 0;
t_xd=text(geometry.X_d_cen(1)+3.0, geometry.X_d_cen(2)+4.0, 'X_d','Color','m','FontSize',15);
t_xd.FontWeight = 'bold';

X_cen_fun = matlabFunction(geometry.poly_X);
fc_X = fcontour(X_cen_fun,'-g', 'LineWidth',2);
fc_X.LevelList = 0;
t_x=text(geometry.X_cen(1)+2.5, geometry.X_cen(2)-2.5, 'X','Color',[0.39,0.83,0.07],'FontSize',15);
t_x.FontWeight = 'bold';

end

