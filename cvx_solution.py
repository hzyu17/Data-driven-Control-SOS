import cvxpy as cvx
import numpy as np

delta = 0.1
nsplits = int(10.0 / delta)
npts = int(nsplits * nsplits)
rho = cvx.Variable([npts, 1], name='rho')
bar_rho = cvx.Variable([npts, 1], name='bar_rho')


def obj_L2():
    x_grid = np.linspace(-5.0, 5.0, nsplits)
    y_grid = np.linspace(-5.0, 5.0, nsplits)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
    X_grid, Y_grid = X_grid.reshape([-1, 1]), Y_grid.reshape([-1,1])
    qx = np.power(X_grid, 4) + np.power(Y_grid, 4)
    qx_rho = sum(qx.T * rho)
    print((bar_rho.T * bar_rho / rho).shape)
    return qx_rho + sum(bar_rho.T * bar_rho / rho)

def construct_constraints():
    rho_grid = cvx.reshape(rho,[nsplits, nsplits])
    bar_rho_grid = cvx.reshape(bar_rho, [nsplits, nsplits])
    
    div_1 = np.diag(np.ones([nsplits]))
    i, j = np.indices(div_1.shape)
    div_1[i==j-1] = -1
    div_1[0, 0] = 0.0

    print(div_1)

if __name__ == '__main__':
    obj_L2()
    construct_constraints()
    x = cvx.Variable()
    y = cvx.Variable()
    
    obj = cvx.Minimize((x-y)**2)

    constraints = [x+y == 1, x-y >= 1]
    prob = cvx.Problem(obj, constraints)
    prob.solve()
    print("status:", prob.status)
    print("optimal value", prob.value)
