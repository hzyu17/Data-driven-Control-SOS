function lc_lqr = generate_local_lqr(dynamics_sym, var_sym, dim_m)
%Generate the matrices associated witht the local lqr controller
%approximation.
lc_lqr.A = double(subs(jacobian(dynamics_sym.F,var_sym.x),var_sym.x,[0;0]));  
lc_lqr.Q = [1 0;0 1]; 
lc_lqr.R = eye(dim_m); 
lc_lqr.N = 0;
lc_lqr.B = dynamics_sym.G;

[lc_lqr.K,lc_lqr.S,lc_lqr.e] = lqr(lc_lqr.A,lc_lqr.B,lc_lqr.Q,lc_lqr.R,lc_lqr.N);

end

