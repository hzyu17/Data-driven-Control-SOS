import nlopt
import cvxpy as cvx
import numpy as np
import matplotlib.pyplot as plt

# grid
delta = 0.5
nsplits = int(10.0 / delta)
npts = nsplits * nsplits
x_grid = np.linspace(-5.0, 5.0, nsplits)
y_grid = np.linspace(-5.0, 5.0, nsplits)
X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
n_constraint = (nsplits-1)*(nsplits-1)

# cvx variables
rho = cvx.Variable([nsplits, nsplits], name='rho')
bar_rho = cvx.Variable([nsplits, nsplits], name='bar_rho')
#rho_grid = cvx.reshape(rho,(nsplits, nsplits))
#bar_rho_grid = cvx.reshape(bar_rho, (nsplits, nsplits))

# initial distribution h_0
h_0 = np.zeros([nsplits, nsplits])
disk = np.multiply(X_grid + 1.5, X_grid + 1.5) + np.multiply(Y_grid - 2.5, Y_grid - 2.5)
disk1 = np.multiply(X_grid - 1.5, X_grid - 1.5) + np.multiply(Y_grid + 3.0, Y_grid + 3.0)
h_0[np.where(disk <= 0.25)] = 1.0
h_0[np.where(disk1<0.25)] = 1.0

# unsafe set condition X_u
bar_lambda = 1
X_u = np.zeros([nsplits, nsplits])
disk = np.multiply(X_grid-0.15, X_grid-0.15) + np.multiply(Y_grid-1.0, Y_grid-1.0)
disk1 = np.multiply(X_grid+0.3, X_grid+0.3) + np.multiply(Y_grid+1.2, Y_grid+1.2)
X_u[np.where(disk <= 0.25)] = 1.0
X_u[np.where(disk1 <= 0.25)] = 1.0
# X_u
plt.contourf(x_grid, y_grid, X_u)
plt.show()
# h_0
plt.contourf(x_grid, y_grid, h_0)
plt.show()

# VanderPol dynamics
f1 = Y_grid
f2 = np.multiply((1.0 - np.multiply(X_grid, X_grid)), Y_grid) - X_grid
plt.quiver(X_grid, Y_grid, f1, f2)
plt.show()
g2 = np.ones(Y_grid.shape)

# state cost
#X_grid_vec, Y_grid_vec = X_grid.reshape([-1, 1]), Y_grid.reshape([-1,1])
qx = np.power(X_grid, 4) + np.power(Y_grid, 4)
qx_rho = cvx.sum(cvx.multiply(qx, rho))

# constraints
def constraints(rho, bar_rho, var_s):
    h_0 = np.ones([nsplits, nsplits])[:-1, :-1]

    rho_f1 = cvx.multiply(rho, f1)
    rho_f2 = cvx.multiply(rho, f2)
    bar_rho_g = cvx.multiply(bar_rho, g2)

    div = np.diag(np.ones([nsplits]))
    i, j = np.indices(div.shape)
    div[i==j-1] = -1
    div_1 = div[:, 1:]
    div_2 = -div[:-1, :]
    
    finite_diff_rhof_x1 = (rho_f1@div_1)[:-1, :] / delta
    finite_diff_rhof_x2 = (div_2@rho_f2)[:, :-1] / delta
    finite_diff_barrho_x2 = (div_2@bar_rho_g)[:, :-1] / delta
    

constraints = []
for i in range(nsplits-1):
    for j in range(nsplits-1):
        constraints.append(finite_diff_rhof_x1[i,j] + finite_diff_rhof_x2[i,j] + finite_diff_barrho_x2[i,j] >= 0)

def l2_cost():
    # L2 cost: state + unsafe + control cost
    cost_Xu = cvx.multiply(bar_lambda, X_u)
    cost_Xu = cvx.sum(cost_Xu)
    #qx_rho# + cost_Xu# + cvx.sum(bar_rho.T * bar_rho / rho)
    l2_cost = qx_rho# + cvx.norm2(cvx.multiply(bar_rho, bar_rho))
    obj = cvx.Minimize(l2_cost)
    return cvx.Problem(obj, constraints)

def l1_cost():
    # l1 cost
    beta = 0.0005
    #l1_cost = qx_rho + cvx.norm1(cvx.multiply(beta, bar_rho))
    var_s = cvx.Variable([nsplits, nsplits])
    l1_cost = cvx.multiply(beta, var_s)
    l1_cost = qx_rho + cvx.sum(l1_cost)
    obj = cvx.Minimize(l1_cost)
    for i in range(nsplits):
        for j in range(nsplits):
            constraints.append(rho[i, j] <= var_s[i, j])
            constraints.append(-rho[i, j] <= var_s[i, j])
            constraints.append(var_s[i, j] >= 0)
    return cvx.Problem(obj, constraints)


if __name__ == '__main__':
    #p = l2_cost()
    p = l1_cost()

    #p.solve(solver=cvx.SCS, verbose=True)
    p.solve(solver=cvx.ECOS_BB, verbose=True)
    
    print("status:", p.status)
    print("optimal value", p.value)
