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
u = sym('u',[0,1],'real');
dx = sym('dx',[2,1],'real');

% model, dxdt = fx + gx*u
fx = [x(2); (1 - x(1)^2) * x(2) - x(1)];
gx = [0; 1];

A = double(subs(jacobian(fx,x),x,[0;0]));
B = gx;
Q = eye(2);
R = 1;
N = 0;
[K_lqr,S,e] = lqr(A,B,Q,R,N);

pvar x1 x2 u
x = [x1; x2];
poly_b = x'*S*x;

% Choice 2. generic bx
% poly_b = x'*x;
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

% coeff_x
c_x = poly2basis(x', Psi);

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

save('polynomials_abc.mat')