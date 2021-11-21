import cvxpy as cvx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys

if len(sys.argv) == 1:
    bar_lambda = 0
else:
    bar_lambda = sys.argv[1]

# grid
delta = 0.2
nsplits = int(10.0 / delta)
npts = nsplits * nsplits
x_grid = np.linspace(-5.0, 5.0, nsplits)
y_grid = np.linspace(-5.0, 5.0, nsplits)
X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
n_constraint = (nsplits-1)*(nsplits-1)

# cvx variables
rho = cvx.Variable([nsplits, nsplits], name='rho')
bar_rho = cvx.Variable([nsplits, nsplits], name='bar_rho')

# initial distribution h_0
h_0 = np.zeros([nsplits, nsplits])
disk = np.multiply(X_grid + 1.5, X_grid + 1.5) + np.multiply(Y_grid - 2.5, Y_grid - 2.5)
disk1 = np.multiply(X_grid - 1.5, X_grid - 1.5) + np.multiply(Y_grid + 3.0, Y_grid + 3.0)

h_0[np.where(disk <= 0.25)] = 1.0
h_0[np.where(disk1 <= 0.25)] = 1.0

# gaussian distribution
sig = 1.0
Sigma = np.eye(2) / sig
MU_1 = np.array([-1.5, 2.5])
MU_2 = np.array([1.5, -3.0])
h_0_gaussian = np.exp(-(np.power(np.multiply(X_grid - MU_1[0], 1/sig), 2) + np.power(np.multiply(Y_grid - MU_1[1], 1/sig), 2))) / (2*np.pi*np.linalg.det(Sigma))
h_0_gaussian += np.exp(-(np.power(np.multiply(X_grid - MU_2[0], 1/sig), 2) + np.power(np.multiply(Y_grid - MU_2[1], 1/sig), 2))) / (2*np.pi*np.linalg.det(Sigma))

# unsafe set condition X_u
X_u = np.zeros([nsplits, nsplits])
disk = np.multiply(X_grid-0.15, X_grid-0.15) + np.multiply(Y_grid-1.0, Y_grid-1.0)
disk1 = np.multiply(X_grid+0.3, X_grid+0.3) + np.multiply(Y_grid+1.2, Y_grid+1.2)
X_u[np.where(disk <= 0.25)] = 1.0
X_u[np.where(disk1 <= 0.25)] = 1.0

# gaussian distribution
sig_Xu = 1.0
Sigma_Xu = np.eye(2) / sig
MU_Xu1 = np.array([0.15, 1.0])
MU_Xu2 = np.array([-0.3, -1.2])
X_u_gaussian = np.exp(-(np.power(np.multiply(X_grid - MU_Xu1[0], 1/sig_Xu), 2) + np.power(np.multiply(Y_grid - MU_Xu1[1], 1/sig_Xu), 2))) / (2*np.pi*np.linalg.det(Sigma_Xu))
X_u_gaussian += np.exp(-(np.power(np.multiply(X_grid - MU_Xu2[0], 1/sig_Xu), 2) + np.power(np.multiply(Y_grid - MU_Xu2[1], 1/sig_Xu), 2))) / (2*np.pi*np.linalg.det(Sigma_Xu))

# plot h_0 + X_u
fig, axs = plt.subplots(2, 1, constrained_layout=True, subplot_kw={"projection": "3d"})

surf_h0 = axs[0].plot_surface(X_grid, Y_grid, h_0_gaussian, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf_Xu = axs[1].plot_surface(X_grid, Y_grid, X_u_gaussian, cmap=cm.coolwarm, linewidth=0, antialiased=False)

axs[0].set_title('h_0_gaussian')
axs[1].set_title('X_u_gaussian')
plt.xlabel("x1")
plt.ylabel("x2")
#plt.show()

# VanderPol dynamics
f1 = Y_grid
f2 = np.multiply((1.0 - np.multiply(X_grid, X_grid)), Y_grid) - X_grid
fig1, axs1 = plt.subplots(1, 1, constrained_layout=True)
axs1.set_title('Van der Pol')
axs1.quiver(X_grid, Y_grid, f1, f2)
#plt.show()
g2 = np.ones(Y_grid.shape)

# state cost qx
qx = np.power(X_grid, 4) + np.power(Y_grid, 4)
fig2, axs2 = plt.subplots(subplot_kw={"projection": "3d"})
axs2.set_title('state cost')
surf = axs2.plot_surface(X_grid, Y_grid, qx, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel("x1")
plt.ylabel("x2")
#plt.show()



rho_f1 = cvx.multiply(rho, f1)
rho_f2 = cvx.multiply(rho, f2)
bar_rho_g1 = cvx.multiply(np.zeros(bar_rho.shape), bar_rho)
bar_rho_g2 = cvx.multiply(bar_rho, g2)

## test divergence by matrix operations                                                                                                                                                   
fx_test = np.array(range(16)).reshape([4,4])
fy_test = np.array(range(1,17)).reshape([4,4])
def test_divergence(fx, fy, i, j):
    matrix_x = np.zeros(fx.shape)
    matrix_x[j,j] = -1
    matrix_x[j+1, j] = 1
    
    construct = np.zeros([fx.shape[0], 1])
    construct[i] = 1

    matrix_y_i = np.zeros(fy.shape)
    matrix_y_i[j,j] = 1
    matrix_y_i[j+1, j] = 1
    
    matrix_y_iplusone = np.zeros(fy.shape)
    matrix_y_iplusone[j,j] = -1
    matrix_y_iplusone[j+1,j] = -1
    
    return (cvx.matmul(cvx.reshape(construct, (fx.shape[0], 1)), cvx.reshape(cvx.matmul(fx[i,:], matrix_x) + cvx.matmul(fx[i+1,:],matrix_x), (1, fx.shape[0]))) + cvx.matmul(cvx.reshape(construct, (fy.shape[0], 1)), cvx.reshape(cvx.matmul(fy[i,:], matrix_y_i) + cvx.matmul(fy[i+1,:],matrix_y_iplusone), (1, fy.shape[0])))) / (2.0 * delta)
    #return (np.outer(construct, np.matmul(fx[i,:], matrix_x) + np.matmul(fx[i+1,:],matrix_x)) + np.outer(construct, np.matmul(fy[i,:], matrix_y_i) + np.matmul(fy[i+1,:],matrix_y_iplusone))) / (2.0 * delta)

#print(test_divergence(fx_test, fx_test, 1, 2))

def div_f(fx, fy):
    ## divergence approx by divergence integration
    div_map = np.zeros([nsplits, nsplits])
    for i in range(1, nsplits-1):
        for j in range(1, nsplits-1):
            div_map = div_map + test_divergence(fx, fy, i, j)
    return div_map

div_frho = div_f(rho_f1, rho_f2)
div_gbar_rho = div_f(bar_rho_g1, bar_rho_g2)
print(div_frho.shape)

# constraints
h_0 = h_0[1:-1, 1:-1]

# constraints on divergence
constraints = []
for i in range(nsplits-2):
    for j in range(nsplits-2):
        constraints.append(div_frho[i,j] + div_gbar_rho[i,j] == h_0_gaussian[i,j])

for i in range(nsplits):
    for j in range(nsplits):
        constraints.append(rho[i, j] >= 0)

dx_squre = np.ones([nsplits, nsplits]) * delta * delta
cost_Xu = cvx.multiply(cvx.multiply(bar_lambda, X_u_gaussian), dx_squre)
cost_Xu = cvx.sum(cost_Xu)
qx_rho = cvx.sum(cvx.multiply(cvx.multiply(qx, rho), dx_squre))

def l2_cost():
    # L2 cost: state + unsafe + control cost
    wx = cvx.Variable([nsplits, nsplits])
    for i in range(nsplits):
        for j in range(nsplits):
            M = cvx.Variable((2, 2), PSD=True)
            constraints.append(M[0,0] == wx[i,j])
            constraints.append(M[0,1] == bar_rho[i,j])
            constraints.append(M[1,0] == bar_rho[i,j])
            constraints.append(M[1,1] == rho[i,j])
    
    #qx_rho# + cost_Xu# + cvx.sum(bar_rho.T * bar_rho / rho)
    l2_cost = qx_rho + cvx.sum(cvx.multiply(wx, dx_squre)) + cost_Xu# + cvx.norm2(cvx.multiply(bar_rho, bar_rho))
    obj = cvx.Minimize(l2_cost)
    return cvx.Problem(obj, constraints)

def l1_cost():
    # l1 cost
    beta = 0.05
    l1_cost = qx_rho + cost_Xu + cvx.norm1(cvx.multiply(cvx.multiply(beta, bar_rho), dx_squre))
    obj = cvx.Minimize(l1_cost)
    return cvx.Problem(obj, constraints)

#def l1_cost():
    # l1 cost
#    beta = 0.0005
#    l1_cost = qx_rho + cvx.norm1(cvx.multiply(beta, bar_rho))
    #var_s = cvx.Variable([nsplits, nsplits], name='var_s')
    #l1_cost = cvx.multiply(beta, var_s)
    #l1_cost = qx_rho + cvx.sum(l1_cost) + cvx.sum(cvx.multiply(rho,cvx.multiply(bar_lambda, X_u)))
#    obj = cvx.Minimize(l1_cost)
    #for i in range(nsplits):
    #    for j in range(nsplits):
    #        constraints.append(bar_rho[i, j] <= var_s[i, j])
    #        constraints.append(-bar_rho[i, j] <= var_s[i, j])
    #        constraints.append(var_s[i, j] >= 0)
#    return cvx.Problem(obj, constraints)

def sys_trj(solved_p):
    for v in solved_p.variables():
        #print(v.name())
        if v.name() == 'rho':
            rho_sol = v.value
        elif v.name() == 'bar_rho':
            bar_rho_sol = v.value
    #print(rho_sol)
    #print(bar_rho_sol)
    
    u_sol = bar_rho_sol / rho_sol
    f2_sol = f2 + np.multiply(g2, u_sol)
    
    fig3, axs3 = plt.subplots(1, 1, constrained_layout=True)
    axs3.set_title('solved system')
    axs3.quiver(X_grid, Y_grid, f1, f2_sol)
    plt.show()


if __name__ == '__main__':
    #p = l2_cost()
    p = l1_cost()

    #p.solve(solver=cvx.SCS, verbose=True)
    #p.solve(solver=cvx.ECOS_BB, verbose=True)
    p.solve(verbose=True)

    print("status:", p.status)
    print("optimal value", p.value)
    
    sys_trj(p)


