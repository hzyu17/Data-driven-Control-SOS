function geo = domain_definition()
% Define the control problem domain
syms x1 x2
% configuration of unsafe sets
num_xu = 2;
num_xinit = 2;
X_u_cen = zeros(num_xu, 2);
X_init_cen = zeros(num_xinit, 2);
r_X = 4.2;
r_Xd = 4.5;
r_xinit = [0.5, 0.5];
r_xu = [0.5, 0.5];
r_xr = 0.4;

%% vdp L2
% vdp: cdc paper
% scene 1
% num_xu = 1;
% num_xinit = 1;
% X_init_cen = [-1.2 2.2];
% X_u_cen(1, :) = [-1.6 0.5];

% scene 2
% X_init_cen = [-1.5 0];
% X_u_cen(1, :) = [1.0, 1.25];

%% VDP L1
% X_init_cen(1, :) = [-1.5 2.5];
% X_init_cen(2, :) = [1.5 -3.0];
% 
% X_u_cen(1, :) = [0.4, 1.0];
% X_u_cen(2, :) = [-0.5, -1.5];

%% VDP L2 2 sets
X_init_cen(1, :) = [1.5 2.5];
X_init_cen(2, :) = [-1.6 -2.5];

X_u_cen(1, :) = [-0.5, 1.2];
X_u_cen(2, :) = [0.5, -1.3];

% 2d linear system
% X_init_cen = [-1.2 2.2];
% X_u_cen(1, :) = [-1.4 1.0];

%% double integrator
% X_init_cen = [-1.2, 1.2];

% X_u_cen = [-1.2, 1.2];
% X_u_cen(1, :) = [0.5 1.2];

%% rantzer's example 2
% X_u_cen(1, :) = [-1.0, -1.0];
% X_init_cen = [1.5, 0];
% % X_init_cen = [-1.5, 1.2];
% r_xu = [0.4];
% r_xinit = 0.5;
% r_xr = 0.5;

%% linear
% X_u_cen(1, :) = [-1.0, -1.0];
% X_init_cen = [-1.5, 1.2];
% r_xu = [0.4];
% r_xinit = 0.5;
% r_xr = 0.5;

%% minus cubic
% X_init_cen = [-1.5 2.5];
% X_u_cen(1, :) = [-1.2 0.5];

%%
X_r_cen = [0 0];
X_d_cen = [0 0];
X_cen = [0 0];

%% 
% domain centers, for OCP integration
geo.X_0.xmin = -r_X;
geo.X_0.xmax = r_X;
geo.X_0.ymin = -r_X;
geo.X_0.ymax = r_X;

geo.X_init = [];
for i_xinit = 1:num_xinit
    X_init.xmin = X_init_cen(i_xinit, 1)-r_xinit(i_xinit);
    X_init.xmax = X_init_cen(i_xinit, 1)+r_xinit(i_xinit);
    X_init.ymin = X_init_cen(i_xinit, 2)-r_xinit(i_xinit);
    X_init.ymax = X_init_cen(i_xinit, 2)+r_xinit(i_xinit);
    geo.X_init = [geo.X_init, X_init];
end

geo.X_excld.xmin = -0.01;
geo.X_excld.xmax = 0.01;
geo.X_excld.ymin = -0.01;
geo.X_excld.ymax = 0.01;

geo.X_r_coord.xmin = X_r_cen(1)-r_xr;
geo.X_r_coord.xmax = X_r_cen(1)+r_xr;
geo.X_r_coord.ymin = X_r_cen(2)-r_xr;
geo.X_r_coord.ymax = X_r_cen(2)+r_xr;

geo.X_u_coord = [];
for i_xu = 1:num_xu
    X_u_coord_i.xmin = X_u_cen(i_xu, 1)-r_xu(i_xu);
    X_u_coord_i.xmax = X_u_cen(i_xu, 1)+r_xu(i_xu);
    X_u_coord_i.ymin = X_u_cen(i_xu, 2)-r_xu(i_xu);
    X_u_coord_i.ymax = X_u_cen(i_xu, 2)+r_xu(i_xu);
    geo.X_u_coord = [geo.X_u_coord, X_u_coord_i];
end

%%
% for plotting
geo.r_xinit = r_xinit;
geo.r_xr = r_xr;
geo.r_xu = r_xu;
geo.r_X = r_X;
geo.r_Xd = r_Xd;

geo.X_init_cen = X_init_cen;
geo.X_u_cen = X_u_cen;
geo.X_r_cen = X_r_cen;
geo.X_d_cen = X_d_cen;
geo.X_cen = X_cen;

%%
% for safety constraints
geo.poly_X = r_X^2 - x1^2 - x2^2;
geo.poly_X_d = r_Xd^2 - x1^2 - x2^2;
geo.poly_X_r = r_xr^2 - (x1-X_r_cen(1))^2 - (x2-X_r_cen(2))^2;

poly_X_init = [];
for i_xinit = 1:num_xinit
    poly_X_init = [poly_X_init, r_xinit(i_xinit)^2 - (x1-X_init_cen(i_xinit, 1))^2 - (x2-X_init_cen(i_xinit, 2))^2];
end
geo.poly_X_init = poly_X_init;

poly_X_u = [];
for i_xu = 1:num_xu
    poly_X_u = [poly_X_u, r_xu(i_xu)^2 - (x1-X_u_cen(i_xu,1))^2 - (x2-X_u_cen(i_xu,2))^2];
end
geo.poly_X_u = poly_X_u;

geo.num_xu = num_xu;
geo.num_xinit = num_xinit;
end

