function geo = domain_definition()
% Define the control problem domain
syms x1 x2

r_xinit = 0.5^2;
r_xr = 0.5^2;
r_xu = 0.5^2;
r_X = 4.0^2;
r_Xd = 4.2^2;

X_init_cen = [1.5 0];
X_u_cen = [-0.55 -1.5];
X_r_cen = [0 0];
X_d_cen = [0 0];
X_cen = [0 0];

% for OCP integration
geo.X_0 = [[-1.0*r_X, r_X, -1*r_X, r_X]];
geo.X_init = [X_init_cen(1)-r_xinit, X_init_cen(1)+r_xinit,  X_init_cen(2)-r_xinit, X_init_cen(2)+r_xinit];
geo.X_excld = [-0.1, 0.1, -0.1, 0.1];
geo.X_r_coord = [X_r_cen(1)-r_xr, X_r_cen(1)+r_xr,  X_r_cen(2)-r_xr, X_r_cen(2)+r_xr];
geo.X_u_coord = [X_u_cen(1)-r_xu, X_u_cen(1)+r_xu,  X_u_cen(2)-r_xu, X_u_cen(2)+r_xu];

% for plotting
geo.r_xinit = r_xinit;
geo.r_xr = r_xr;
geo.r_xu = r_xu;
geo.r_X = 5.0^2;
geo.r_Xd = 5.2^2;

geo.X_init_cen = X_init_cen;
geo.X_u_cen = X_u_cen;
geo.X_r_cen = X_r_cen;
geo.X_d_cen = X_d_cen;
geo.X_cen = X_cen;

% for safety constraints
geo.X = r_X - x1^2 - x2^2;
geo.X_d = r_Xd - x1^2 - x2^2;
geo.X_init = geo.r_xinit - (x1-X_init_cen(1))^2 - (x2-X_init_cen(2))^2;
geo.X_r = geo.r_xr - (x1-X_r_cen(1))^2 - (x2-X_r_cen(2))^2;
geo.X_u = geo.r_xu - (x1-X_u_cen(1))^2 - (x2-X_u_cen(2))^2;

end

