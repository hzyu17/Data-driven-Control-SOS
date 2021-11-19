%% CVX solve for the optimal control problem
% L-S problem to verigy that cvx runs normally
m=16;n=8;
A=randn(m,n);
b=randn(m,1);

cvx_begin
    variable x(n)
    minimize ( norm(A*x-b))
cvx_end
