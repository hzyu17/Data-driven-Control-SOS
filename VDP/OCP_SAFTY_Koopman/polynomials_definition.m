function polynomials_definition(coeff_q)
% generating polynomials

% VDP
total_deg = 14;
pvar x1 x2
x = [x1; x2];
nx = length(x);

nPsi = nchoosek(total_deg+nx, total_deg);
Psi = monomials([x1 x2], [0:total_deg]);
Psi_sym = p2s(Psi);

% optimal control polynomials
deg_a = 4;

% dimension of the control input c(x)
dim_m = 2;
deg_c = [6, 6];

if length(deg_c) ~= dim_m
    disp('wrong dimension of inputs')
    keyboard
end

max_degc = max(deg_c);
deg_w = max_degc*2 - deg_a;
Alph = 4; 

% define for ax cx wx
Qa = nchoosek(deg_a+nx, nx);
c_a = mpvar('c_a', [Qa, 1]);
c_a = [c_a; zeros(nPsi - Qa, 1)];
c_a_sym = p2s(c_a);
poly_a = c_a'* Psi;
poly_a_sym = c_a_sym' * Psi_sym;

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
    c_c_sym = [c_c_sym; {c_c_i_sym}];
    poly_c_i = c_c_i' * Psi;
    poly_c = [poly_c; poly_c_i];
    poly_c_sym_i = c_c_i_sym' * Psi_sym;
    poly_c_sym = [poly_c_sym; poly_c_sym_i];
end

Qw = nchoosek(deg_w+nx, nx);
c_w = mpvar('c_w', [Qw, 1]);
c_w = [c_w; zeros(nPsi - Qw, 1)];
c_w_sym = p2s(c_w);
poly_w_sym = c_w_sym' * Psi_sym;

pvar x1 x2
x = [x1; x2];

[fx, gx] = dynamics_vdp();
dim_m = size(gx);
dim_m = dim_m(2);

% local controller design using LQR
syms x1 x2
x = [x1;x2];
A = double(subs(jacobian(fx,x),x,[0;0]));  
Q = [1 0;0 1]; 
R = eye(dim_m); 
N = 0;
B = gx;

[K,S,e] = lqr(A,B,Q,R,N);
I = [1 0;0 1];
if strcmp(class(x1), 'polynomial') == 1
    poly_b = x.'*S*x;
    poly_b_sym = p2s(poly_b);
else
    poly_b_sym = x.'*S*x;
end

% define for qx
syms x1 x2
% poly_q_sym = 3.5*x1^2 + x2^2;
%inverted pendulum
poly_q_sym = x1^2 + x2^2;

% define for abx, acx
poly_bc_sym = [];
poly_ab_sym = poly_a_sym * poly_b_sym;
for i_c = 1:dim_m
    poly_bc_sym = [poly_bc_sym; poly_b_sym * poly_c_sym(i_c)];
end

%%
syms x1 x2
x = [x1; x2];
c_bc_sym = [];
for i_c = 1:dim_m
    c_bc_sym = [c_bc_sym; getSymCoeffs(poly_bc_sym(i_c), x, Psi)];
end
c_ab_sym = getSymCoeffs(poly_ab_sym, x, Psi);

save('polynomials_def')

end

