deg_a  = 2;
deg_b = 2;
deg_c = 4;
deg_total = 8;
pvar x1 x2
x = [x1;x2];
bx = x' * x;

opt.round = 1e-6;
opt.errnrm = 1e-8;

[C_a_Psi, C_cx_Psi, C_ab_Psi, C_bc_Psi, Psi] = gen_polynomial(x, bx, deg_a, deg_c, deg_total, opt);

load('external/VdP gEDMD/gEDMD/result_gEDMD.mat');