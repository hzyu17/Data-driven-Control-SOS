%% Notes and logs
% 01/12/2021 cellfun, matlabFunction is much more faster than subs() for
% symbolic substitutions


close all 
clear all

% ------------------------------
% Variables and dimensions 
% ------------------------------
pvar x1 x2 u
x = [x1;x2];
nx = length(x); % variable dimension
m = length(u); % input dimension

%% 
% -------------------------
% Polynomial degrees
% -------------------------
deg_c = 4;
deg_a = 2;
% highest_order = max(deg_c+2, deg_a+deg_c)
highest_order = 8
nPsi = nchoosek(highest_order+nx, nx);
alpha = 4;

% -------------------------------------------
% Psi: common basis for all polynomials
% -------------------------------------------
Psi = monomials(x, 0:highest_order);
if size(Psi, 1) ~= nPsi
    fprintf('wrong number of monomials Psi!');
end

% ---------------------------------------
% Setting bx as a known polynomial
% ---------------------------------------

% Choice 1. local controller design using LQR, CLF for bx
% variable definition
x = sym('x',[2,1],'real');
u = sym('u',[1,1],'real');
model_dyn = sym('dx',[2,1],'real');
t = sym('t', 0, 'real'); % never used

% model, dxdt = fx + gx*u
fx = model_dynamics_VDP(t, x, u);
gx = [0; 1];

A = double(subs(jacobian(fx,x),x,[0;0]));
B = gx;
Q = eye(2);
R = 1;
N = 0;
[K_lqr,S,e] = lqr(A,B,Q,R,N);

pvar x1 x2 u
x = [x1; x2];
% poly_b = x'*S*x;

% Choice 2. generic bx
poly_b = x'*x;
c_b = poly2basis(poly_b, Psi);

% ----------------------------------------------
% c(x) and a(x) is the unknown polynomial
% ----------------------------------------------

% cx
Qc = nchoosek(deg_c + nx, nx);
c_c = mpvar('c', [Qc, 1]);
% c_c(1) = 0;
z = zeros(nPsi-Qc, 1);
p = polynomial(z);
c_c = [c_c;p];
poly_c = c_c' * Psi;

% ax
Qa = nchoosek(deg_a + nx, nx);
c_a = mpvar('a', [Qa, 1]);
c_a = [c_a; zeros(nPsi-Qa,1)];
poly_a = c_a' * Psi;

% --------------------------------------------------------
% Get mixed polynomial coefficients using symbolic
% --------------------------------------------------------
syms x1 x2
x = [x1;x2];

% 1. Coefficients of  abx
poly_ab_sym = p2s(poly_a * poly_b);
c_ab = getSymCoeffs(poly_ab_sym, x, Psi);

% 2. Coefficients of bcx
poly_bc_sym = p2s(poly_c * poly_b);
c_bc = getSymCoeffs(poly_bc_sym, x, Psi);

% 3. Coefficients of x
pvar x1 x2 u
x = [x1;x2];
c_x = poly2basis(x', Psi);

%%
% ========================
% sampling gEDMD for X's and X_dot's
% ========================

% generate random initial points and collector for the next point

u = 0; % input
t = 0; % created but not used 

for i = 1 : m+1
    Xinit = -5.0 + 10.0 * rand(nx, 10000); % random initial points
    
    if i == 2
        u = 1;
    end
    
    % collecting for the u_i simulation
    Xinit_i = [];
    Yinit_i = [];
    % solve ODE and sample
    if strcmp(class(x), 'polynomial')
        x1 = p2s(x1);
        x2 = p2s(x2);
        x = [x1; x2];
    end
    
    %% gEDMD
    
    % Generate basis Psi and its time derivative expression
    Psi_sym = p2s(Psi);
    grad_Psi = jacobian(Psi_sym);
    model_dyn = model_dynamics_VDP(t, x, u);
    Psi_dt = grad_Psi * model_dyn;    
    
    % Get initial data Psi(x0)
    data_cell = mat2cell(Xinit, size(Xinit, 1), ones(1, size(Xinit, 2)));
    
    MFun_Psi = matlabFunction(Psi_sym, 'Vars', {[x1;x2]});
    Xi = cell2mat(cellfun(MFun_Psi, data_cell, 'UniformOutput', 0));
    
    % Get derivative data on the basis dPsi(x0)
    
    % Important! Use matlabFunction and cellfun to evaluate large
    % scales of data is much more efficient than subs!
    MFun_Psi_dt = matlabFunction(Psi_dt, 'Vars', {[x1;x2]});
    dXi =  cell2mat(cellfun(MFun_Psi_dt, data_cell, 'UniformOutput', 0));
    
    % approximation
    M_i = size(Xi, 2);
    A = 1/M_i * dXi * Xi';
    B = 1/M_i * Xi * Xi';
    Li = A * pinv(B);
    figure;
    spy(Li)
    title('Sparsity visualization of the approximated infinitesimal matrix L')
    
    if i==1
        X = Xi;
        Y = dXi;
        L = Li;
    else
        X = {X,Xi};
        Y = {Y,dXi};
        L = {L, Li};
    end

end
%%
formatOut = 'mm_dd_yy_hh_MM_SS';
ch_date = datestr(now,formatOut);
ch_file = ['data/gEDMD_res_', ch_date,'.mat'];
save(ch_file)