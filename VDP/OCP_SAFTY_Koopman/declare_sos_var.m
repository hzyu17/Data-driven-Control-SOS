function [sos_program, poly_sym] = declare_sos_var(sos_program, dim_x, var_name, Psi, degree, sparse)
nPsi = length(Psi);
Q = nchoosek(degree+dim_x, dim_x);
c_poly = mpvar(var_name, [Q, 1]);
c_poly = [c_poly(1:Q); zeros(nPsi-Q,1)];
poly_sym = p2s(transpose(c_poly) * Psi);
c_poly_sym = p2s(c_poly);
sos_program = sosdecvar(sos_program, c_poly_sym(1:Q));
if strcmp(sparse, 'true')
    sos_program = sosineq(sos_program, poly_sym,'sparse');
else
    sos_program = sosineq(sos_program, poly_sym);
end
end

