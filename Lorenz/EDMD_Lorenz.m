% EDMD method to approximate Koopman operator
% ------------------
% Main content:
% 1. Generate all polynomials ax, bx, cx, abx, bcx, and basis Psi;
% 2. Calculate all matrices for the Koopman: K, L, P.
% ------------------
% HISTORY:
% 12/01 wrap polynomial basis coefficient
% 12/02 implement polynomial coefficient generator function
% 12/02 test the result with the new method generating polynomials: same as
% the symbolic way. Symbolic way is easier.
% 12/03 Read sampling data from Hyungjin's data, and do the EDMD, to test
%           my sampling methods -- result: no the sampling's fault.
% 12/03 Test with matlab symbolic and the MatlabFunction for Psi, to take
%           values. To see if it is more accurate than polynomial subs().
%           Result: greatly changed the vector field of VdP! towards the
%           correct one except for one divergence.
% 12/03 Test the one step inverse, instead of sum all the Psi_X*Psi_X' and
%           Psi_Y*Psi_X'. Result: not changing much.
% 12/04 Test the normalization of the data points in terms of time. To see
%           the impact: sosopt run with warnings and cannot proceed.

close all
clear all
clc
% ------------------------
% Load sampling data 
% ------------------------
load('sampling_data.mat')

%% Test with Hyungjin's sample data
% load('external/VdP gEDMD/sampling/samples.mat');
% X{1} = Xe1;
% X{2} = Xe2;
% Y{1} = Ye1;
% Y{2} = Ye2;
%%
% ------------------------------
% Variables and dimensions 
% ------------------------------
pvar x1 x2 x3 u
x = [x1;x2;x3];
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
x = sym('x',[3,1],'real');
u = sym('u',[1,1],'real');
dx = sym('dx',[2,1],'real');

% model, dxdt = fx + gx*u
rho_Lor = 28;
sigma_Lor = 10;
beta_Lor = 8/3;
fx = [sigma_Lor * (x(2)-x(1)); 
    x(1) * (rho_Lor - x(3)) - x(2) + u;
    x(1) * x(2) - beta_Lor * x(3)
    ];
gx = [0; 1; 0];

A = double(subs(jacobian(fx,x),x,[0;0;0]));
B = gx;
Q = eye(nx);
R = 1;
N = 0;
[K_lqr,S,e] = lqr(A,B,Q,R,N);

pvar x1 x2 x3 u
x = [x1; x2; x3];
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
syms x1 x2 x3
x = [x1;x2;x3];

% 1. Coefficients of  abx
poly_ab_sym = p2s(poly_a * poly_b);
c_ab = getSymCoeffs(poly_ab_sym, x, Psi);

% 2. Coefficients of bcx
poly_bc_sym = p2s(poly_c * poly_b);
c_bc = getSymCoeffs(poly_bc_sym, x, Psi);

% 3. Coefficients of x
pvar x1 x2 x3 u
x = [x1;x2;x3];
c_x = poly2basis(x', Psi);

% %% 
% % -------------------------------------------------------------
% % Test polynomial coefficient generator function, 
% % and compare with the symbolic method.
% % -------------------------------------------------------------
% % options: round-off threshold, and error norm threshold
% opt.round = 1e-6;
% opt.errnrm = 1e-6;
% % call the function
% [C_a, C_c, C_ab, C_bc, Psi] = gen_polynomial(x, poly_b, deg_a, deg_c, highest_order, opt);
% % testing: comparason. Result: all the coefficients are the same
% if sum(p2s(C_a - c_a)) ~= 0 || sum(p2s(C_c - c_c))~=0 || sum(p2s(C_ab - c_ab))~=0 || sum(p2s(C_bc - c_bc)) ~= 0 || sum(p2s(Psi - Phix)) ~=0
%     fprintf('There are differences between the two methods! Please check the polynomial coefficient generation mthod.')
% else
%     fprintf('The coefficients got from the two methods: symbolic transition, and polynomial transition, are the same!')
% end

%%
% ----------
% EDMD
% ----------
% A and B matrices
A = cell(m+1, 1);
B = cell(m+1, 1);
Xis = {};
Yis = {};

%%
% --------------------------------------------------------------
% Test with symbolic representation and MatlabFunction
% --------------------------------------------------------------
x = sym('x', [3,1], 'real');
Psi_sym = sym({});
for i1 = 1: length(Psi)
    Psi_sym(i1,1) = p2s(Psi(i1));
end
% The Psi function. Will apply this function to all the data points.
Psi_fun = matlabFunction(Psi_sym, 'Vars', {x});
for i1 = 1:m+1
    Ti = size(X{i1},2);
    X_cell = mat2cell(X{i1}, size(X{i1},1), ones(Ti, 1));
    Y_cell = mat2cell(Y{i1}, size(Y{i1},1), ones(Ti, 1));
    
    Psi_Xi{i1} = cell2mat( cellfun(Psi_fun, X_cell, 'UniformOutput', false) );
    Psi_Yi{i1} = cell2mat( cellfun(Psi_fun, Y_cell, 'UniformOutput', false) );
    
    %%
    % ----------------------------------------------
    % Sum all the Psi_X*Psi_X' and Psi_Y*Psi_X'
    % ----------------------------------------------
    XX = zeros(Psi.matdim(1), Psi.matdim(1));
    YX = zeros(Psi.matdim(1), Psi.matdim(1));
    
    for i2 = 1 : Ti
        XX = XX + Psi_Xi{i1}(:, i2) * Psi_Xi{i1}(:, i2)';
        YX = YX + Psi_Yi{i1}(:, i2) * Psi_Xi{i1}(:, i2)';
    end
    
    %%
    % --------------------------------------------------
    % Inverte directly Psi_X*Psi_X' and Psi_Y*Psi_X'
    % --------------------------------------------------
%     XX = Psi_Xi{i1} * Psi_Xi{i1}';
%     YX = Psi_Yi{i1} * Psi_Xi{i1}';
    %%
    % ----------------------------------------
    % Test K matrix after l1-normalization
    % ----------------------------------------
%     norm_PsiX = vecnorm(Psi_Xi{i1}, 1, 2);
%     norm_PsiY = vecnorm(Psi_Yi{i1}, 1, 2);
%     
%     Psi_X_n = Psi_Xi{i1} ./ norm_PsiX;
%     Psi_Y_n = Psi_Yi{i1} ./ norm_PsiY;
%     
%     XX = Psi_X_n * Psi_X_n';
%     YX = Psi_Y_n * Psi_X_n';
    %%
%     if length(A) <= 0
    A{i1} = XX ./ Ti;
    B{i1} = YX ./ Ti;
%     else
%         A =  {A, XX ./ Ti};
%         B =  {B, YX ./ Ti};
%     end

end

% %% 
% % --------------------------------------------------
% % Using polynomial subs to get Psi_X and Psi_Y
% % --------------------------------------------------
% for i = 1:m+1
%     Xi = X{i};
%     Yi = Y{i};
%     Ti = length(Xi);
%     
%     XX = zeros(Psi.matdim(1), Psi.matdim(1));
%     XY = zeros(Psi.matdim(1), Psi.matdim(1));
% 
%     for iter = 1:Ti
%         x_step = Xi(1:2, iter);
%         y_step = Yi(1:2, iter);
%         px_step = subs(Psi, x, x_step);
%         px_step = px_step.coefficient';
%         py_step = subs(Psi, x, y_step);
%         py_step = py_step.coefficient';
%         XX = XX + px_step * px_step';
%         XY = XY + py_step * px_step';
% %         XY = XY + px_step * py_step';
%     end
%     
%     A{i} = XX ./ Ti;
%     B{i} = XY ./ Ti;
% end

%% 
% ----------------------------------------------
% Ki matrices: inifitesimal approximator
% ----------------------------------------------
clc
K = cell(m+1, 1);
for i1 = 1:m+1
    A_pinv = pinv(A{i1});
    K_this = B{i1} * A_pinv;
    %% if l1-normalized K
%     K_this = K_this .* norm_PsiY;
%     K_this = K_this ./ norm_PsiX;
    K{i1} = K_this;
    % Draw sparsity map
    figure;
    spy(K_this);
    title('Sparsity map of the estimated Koopman generator');
%         K = {K, A_pinv * B{i}};
end

%% 
% --------------------------------------------
% L0, L1 for the dynamics approximator
% --------------------------------------------
clc
% L
for i=1:m+1
    if i == 1
        L = (K{1} - eye(nPsi)) ./ dt;
    else
        L = {L, (K{i} - K{1}) ./ dt;};
    end
end

%% 
% ------------------------------------------------------
% Testing: correct value of F and G using L0 and L1
% ------------------------------------------------------
% c_x' * L{1} * Phix: [x2; (1-x1^2)x2-x1]
% c_x' * L{2} * Phix: [0; 1]

pvar x1 x2 x3

%%
% -----------------------------------------------
% P matrices for divergence approximator
% -----------------------------------------------
for j = 1:m+1
    if j == 1
        L_ = L{j};
    else
%         L_ = L{j} - L{j-1};
        L_ = L{j};
    end
    poly_div = c_x' * L_ * Psi;

    % j == 1: div(F) = div([K0(x1), ..., K0(xn)]) = div(c_x' * L0 * Phix)
    % j >=2: div(Gj) = div(c_x' * Lj * Phix)
    div = diff(poly_div(1), x1) + diff(poly_div(2), x2) + diff(poly_div(3), x3);
    
    % polynomial coefficients of this divergence term
    % c_div' * Phix is the div(F) or div(Gj) for [x1, x2]
    
    % j=1: 1-x1^2
    % j=2: 0
    div = clean_polynomial(div, 1e-9);
    
    if j == 1
        P = L{1} + div .* eye(nPsi);
    else
        P = {P, L{j} + div .* eye(nPsi)};
    end
end
%%
% --------------------
% Print out for info
% --------------------
poly_a
poly_c
poly_b
highest_order
%%
% ------------
% Save data
% ------------
save('EDMD_res.mat')