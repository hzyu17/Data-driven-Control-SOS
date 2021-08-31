function polynomials_definition(coeff_q)
% generating polynomials

total_deg = 16;
pvar x1 x2
x = [x1; x2];
nx = length(x);
nPsi = nchoosek(total_deg+nx, total_deg);
Psi = monomials([x1 x2], [0:total_deg]);
Psi_sym = p2s(Psi);

% optimal control polynomials
deg_a = 4;
deg_c = 6;
deg_w = deg_c*2 - deg_a;
Alph = 4; 

% define for ax cx wx
Qa = nchoosek(deg_a+nx, nx);
c_a = mpvar('c_a', [Qa, 1]);
c_a = [c_a; zeros(nPsi - Qa, 1)];
c_a_sym = p2s(c_a);
poly_a = c_a'* Psi;
poly_a_sym = c_a_sym' * Psi_sym;

Qc = nchoosek(deg_c+nx, nx);
c_c = mpvar('c_c', [Qc, 1]);
c_c = [c_c; zeros(nPsi - Qc, 1)];
c_c_sym = p2s(c_c);
poly_c = c_c' * Psi;
poly_c_sym = c_c_sym' * Psi_sym;

% Qc1 = nchoosek(deg_c+nx, nx);
% c_c1 = mpvar('c_c1', [Qc1, 1]);
% c_c1 = [c_c1; zeros(nPsi - Qc1, 1)];
% c_c1_sym = p2s(c_c1);
% poly_c1 = c_c1' * Psi;
% poly_c1_sym = c_c1_sym' * Psi_sym;
% 
% c_cs = [c_c; c_c1];
% c_c_syms = [c_c_sym; poly_c1_sym];
% poly_cs = [poly_c; poly_c1];
% poly_c_syms = [poly_c_sym; poly_c1_sym];

Qw = nchoosek(deg_w+nx, nx);
c_w = mpvar('c_w', [Qw, 1]);
c_w = [c_w; zeros(nPsi - Qw, 1)];
c_w_sym = p2s(c_w);
poly_w_sym = c_w_sym' * Psi_sym;

pvar x1 x2
x = [x1; x2];

[fx, gx] = dynamics_vdp();

% local controller design using LQR
A = double(subs(jacobian(fx,x),x,[0;0]));  Q = [1 0;0 1]; R = 1; N = 0;
B = gx;

[K,S,e] = lqr(A,B,Q,R,N);
I = [1 0;0 1];
poly_b = x'*S*x;
poly_b_sym = p2s(poly_b);

% define for qx
syms x1 x2
poly_q_sym = 3.5*x1^2 + x2^2;

% define for abx, acx
poly_ab_sym = poly_a_sym * poly_b_sym;
poly_bc_sym = poly_b_sym * poly_c_sym;

%%
syms x1 x2
x = [x1; x2];
c_bc_sym = getSymCoeffs(poly_bc_sym, x, Psi);
c_ab_sym = getSymCoeffs(poly_ab_sym, x, Psi);

save('polynomials_def')

end

