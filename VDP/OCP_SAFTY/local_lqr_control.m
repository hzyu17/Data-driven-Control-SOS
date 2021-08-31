function F = local_lqr_control(local_lqr, pt, x)
% The locally linearized lqr controller
A = local_lqr.A;
B = local_lqr.B;
K = local_lqr.K;

F = double(subs(A * x - B*K*x, x, pt));
aa = 1;

end

