function polynomials_definition(option)
% generating polynomials
dynamics_option = option.dynamics_option;

total_deg = option.degrees.total_deg;
deg_a = option.degrees.deg_a;
deg_c = option.degrees.deg_c;
Alph = option.degrees.Alph;
dim_m = option.degrees.dim_m;
deg_q = option.degrees.deg_q;
poly_b_type = option.poly_b_type;

pvar x1 x2
var_poly.x = [x1; x2];
var_poly.x1 = x1;
var_poly.x2 = x2;
nx = length(var_poly.x);

nPsi = nchoosek(total_deg+nx, total_deg);
Psi = monomials(var_poly.x, [0:total_deg]);
Psi_sym = p2s(Psi);

if length(deg_c) ~= dim_m
    disp('wrong dimension of inputs')
    keyboard
end

max_degc = max(deg_c);
deg_w = max_degc*2 - deg_a;

% define for ax cx wx
Qa = nchoosek(deg_a+nx, nx);
c_a = mpvar('c_a', [Qa, 1]);
c_a = [c_a; zeros(nPsi - Qa, 1)];
c_a_sym = p2s(c_a);
% set the constant term in ax is 1: increase feasratio
% c_a_1 = sym(1);
% c_a_sym = [c_a_1; c_a_sym(2:end)];

poly_a = transpose(c_a)* Psi;
poly_a_sym = transpose(c_a_sym) * Psi_sym;

Qc = [];
c_c = [];
c_c_sym = [];
poly_c = [];
poly_c_sym = [];

for i_c = 1:dim_m
    Qc_i = nchoosek(deg_c(i_c)+nx, nx);
    Qc = [Qc, Qc_i];
    var_name = ['c_c', num2str(i_c)];
    c_c_i = mpvar(var_name, [Qc_i, 1]);
    c_c_i = [c_c_i; zeros(nPsi - Qc_i, 1)];
    c_c = [c_c, {c_c_i}];
    c_c_i_sym = p2s(c_c_i);
    % set the constant term in cx is 1: not working feasratio = -1.0097
%     c_c_i_1 = sym(1);
%     c_c_i_sym = [c_c_i_1; c_c_i_sym(2:end)];
    c_c_sym = [c_c_sym; {c_c_i_sym}];
    poly_c_i = transpose(c_c_i) * Psi;
    poly_c = [poly_c; poly_c_i];
    poly_c_sym_i = transpose(c_c_i_sym) * Psi_sym;
    poly_c_sym = [poly_c_sym; poly_c_sym_i];
end

Qw = nchoosek(deg_w+nx, nx);
c_w = mpvar('c_w', [Qw, 1]);
c_w = [c_w; zeros(nPsi - Qw, 1)];
c_w_sym = p2s(c_w);

poly_w_sym = transpose(c_w_sym) * Psi_sym;

dynamics_sym = dynamics_definition('sym', dynamics_option);

dim_m = size(dynamics_sym.G);
dim_m = dim_m(2);

% local controller design using LQR
syms x1 x2
var_sym.x = [x1;x2];
var_sym.x1 = x1;
var_sym.x2 = x2;

% local lqr
if strcmp(poly_b_type, 'lqr')
    lc_lqr = generate_local_lqr(dynamics_sym, var_sym, dim_m);
    poly_b = var_poly.x'*lc_lqr.S*var_poly.x;
    poly_b_sym = var_sym.x.'*lc_lqr.S*var_sym.x;
elseif strcmp(poly_b_type, 'quadratic')
    poly_b = var_poly.x'*var_poly.x;
    poly_b_sym = var_sym.x.'*var_sym.x;
end

% define for qx
% poly_q_sym = 3.5*x1^2 + x2^2;
%inverted pendulum
% poly_q_sym = var_sym.x1^2 + var_sym.x2^2;

% L1 data driven
poly_q_sym = var_sym.x1^deg_q + var_sym.x2^deg_q;

% define for abx, acx
poly_bc_sym = [];
poly_ab_sym = poly_a_sym * poly_b_sym;
for i_c = 1:dim_m
    poly_bc_sym = [poly_bc_sym; poly_b_sym * poly_c_sym(i_c)];
end

%%
c_bc_sym = cell(1, dim_m);
c_bc = cell(1, dim_m);
for i_c = 1:dim_m
    c_bc_i = getSymCoeffs(poly_bc_sym(i_c), var_sym.x, Psi);
    c_bc{1, i_c} = c_bc_i;
    c_bc_sym_i = p2s(c_bc_i);
    c_bc_sym{1,i_c} = c_bc_sym_i;
end

c_ab = getSymCoeffs(poly_ab_sym, var_sym.x, Psi);
c_ab_sym = p2s(c_ab);
% poly_ab = poly_a * poly_b;
% c_ab = poly2basis(poly_ab, Psi);

assume(c_a_sym, 'real');
for i_c = 1:dim_m
    c_c_i = c_c_sym(i_c);
    assume(c_c_i{1:end}, 'real');
end
assume(c_w_sym, 'real');

save('polynomials_def')

end

