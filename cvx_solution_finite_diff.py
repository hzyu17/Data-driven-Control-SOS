import cvxpy as cvx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys

if len(sys.argv) == 1:
    bar_lambda = 0
else:
    bar_lambda = sys.argv[1]

## grid
delta = 0.1
nsplits = int(10.0 / delta)
npts = nsplits * nsplits
x_grid = np.linspace(-5.0, 5.0, nsplits)
y_grid = np.linspace(-5.0, 5.0, nsplits)
X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
n_constraint = (nsplits-1)*(nsplits-1)
dx_square = np.ones([nsplits, nsplits]) * delta * delta

## cvx variables
rho = cvx.Variable(shape=(nsplits, nsplits), name='rho', nonneg=True)
bar_rho = cvx.Variable(shape=(nsplits, nsplits), name='bar_rho')

## initial distribution h_0
# center: [-1.5, 2.5], [1.5, -3.0]
# indicator function
h_0 = np.zeros([nsplits, nsplits])
disk = np.multiply(X_grid + 1.5, X_grid + 1.5) + np.multiply(Y_grid - 2.5, Y_grid - 2.5)
disk1 = np.multiply(X_grid - 1.5, X_grid - 1.5) + np.multiply(Y_grid + 3.0, Y_grid + 3.0)

h_0[np.where(disk <= 0.25)] = 1.0
h_0[np.where(disk1 <= 0.25)] = 1.0

# gaussian distribution
sig = 1.0
Sigma = np.eye(2) / sig
MU_0 = np.array([0.0, 0.0])
MU_1 = np.array([-1.5, 2.5])
MU_2 = np.array([1.5, -3.0])
h_0_gaussian = np.exp(-(np.power(np.multiply(X_grid - MU_1[0], 1/sig), 2) + np.power(np.multiply(Y_grid - MU_1[1], 1/sig), 2))) / (2*np.pi*np.linalg.det(Sigma))
#h_0_gaussian += np.exp(-(np.power(np.multiply(X_grid - MU_2[0], 1/sig), 2) + np.power(np.multiply(Y_grid - MU_2[1], 1/sig), 2))) / (2*np.pi*np.linalg.det(Sigma))

## unsafe set condition X_u
# center: [0.15, 1.0], [-0.3, -1.2]
# indicator function
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
axs[1].set_title('X_u')
plt.xlabel("x1")
plt.ylabel("x2")

## VanderPol dynamics
f1 = Y_grid
f2 = np.multiply((1.0 - np.multiply(X_grid, X_grid)), Y_grid) - X_grid
fig1, axs1 = plt.subplots(1, 1, constrained_layout=True)
axs1.set_title('Van der Pol')
axs1.quiver(X_grid, Y_grid, f1, f2)
g2 = np.ones(Y_grid.shape)

## state cost qx
#X_grid_vec, Y_grid_vec = X_grid.reshape([-1, 1]), Y_grid.reshape([-1,1])
qx = np.power(X_grid, 2) + np.power(Y_grid, 2)
#fig2, axs2 = plt.subplots(subplot_kw={"projection": "3d"})
#axs2.set_title('state cost')
#surf = axs2.plot_surface(X_grid, Y_grid, qx, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#plt.xlabel("x1")
#plt.ylabel("x2")

## constraints
h_0 = h_0[:-1, :-1]

rho_f1 = cvx.multiply(rho, f1)
rho_f2 = cvx.multiply(rho, f2)
#bar_rho_g1 = cvx.multiply(np.zeros(bar_rho.shape), bar_rho)
bar_rho_g2 = cvx.multiply(bar_rho, g2)

## divergence approx by finite difference
div = np.diag(np.ones([nsplits]))
i, j = np.indices(div.shape)
div[i==j-1] = -1
div_1 = div[:, 1:]
#print("div1")
#print(div_1)
div_2 = -div[:-1, :]
#print("div2")
#print(div_2)

numerator = np.ones([nsplits-1, nsplits-1]) / delta

finite_diff_rhof_x1 = cvx.multiply((rho_f1@div_1)[:-1, :], numerator)
finite_diff_rhof_x2 = cvx.multiply((div_2@rho_f2)[:, :-1], numerator)
finite_diff_barrho_x2 = cvx.multiply((div_2@bar_rho_g2)[:, :-1], numerator)

h_0_gaussian = h_0_gaussian[:nsplits-1, :nsplits-1]
# constraints on divergence
constraints = []
constraints.append(finite_diff_rhof_x1 + finite_diff_rhof_x2 + finite_diff_barrho_x2 == h_0_gaussian)
#for i in range(nsplits-1):
#    for j in range(nsplits-1):
#        constraints.append(finite_diff_rhof_x1[i,j] + finite_diff_rhof_x2[i,j] + finite_diff_barrho_x2[i,j] == h_0_gaussian[i,j])
        #constraints.append(bar_rho[i,j] <= 0.5*rho[i,j])
        #constraints.append(finite_diff_rhof_x1[i,j] + finite_diff_rhof_x2[i,j] + finite_diff_barrho_x2[i,j] >= 0)
        #print(div_frho[i,j])
        #print(div_gbar_rho[i,j])
        #constraints.append(div_frho[i,j] + div_gbar_rho[i,j] == h_0[i,j])
#constraints.append(rho >=0)

#for i in range(nsplits):
#    for j in range(nsplits):
#        constraints.append(rho[i, j] >= 0)

#cost_qx_rho = cvx.norm(cvx.multiply(cvx.multiply(qx, rho), dx_square), 1)
cost_qx_rho = cvx.sum(cvx.multiply(cvx.multiply(qx, rho), dx_square))
cost_Xu = cvx.sum(cvx.multiply(cvx.multiply(bar_lambda, X_u_gaussian), dx_square))

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
    l2_cost = cost_qx_rho + cvx.sum(cvx.multiply(wx, dx_square)) + cost_Xu # + cvx.norm2(cvx.multiply(bar_rho, bar_rho))
    obj = cvx.Minimize(l2_cost)
    return cvx.Problem(obj, constraints)

def l1_cost():
    # l1 cost
    beta = 0.0005
    l1_cost = cost_qx_rho + cost_Xu + cvx.norm(cvx.multiply(cvx.multiply(beta, bar_rho), dx_square),1)
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

## solution and solved trajectory
def sys_trj(solved_p):
    for v in solved_p.variables():
        #print(v.name())
        if v.name() == 'rho':
            rho_sol = v.value
        elif v.name() == 'bar_rho':
            bar_rho_sol = v.value
    #print(rho_sol)
    #print(bar_rho_sol)
    print("bar rho")
    print(bar_rho_sol)
    print("rho")
    print(rho_sol)
    u_sol = bar_rho_sol / rho_sol
    print("u_sol")
    print(u_sol)
    f2_sol = f2 + u_sol
    
    ## visualize the solved div(rho_f + bar_rho_g)
    rho_f1_sol = np.multiply(rho_sol, f1)
    rho_f2_sol = np.multiply(rho_sol, f2)
    barrho_g2_sol = np.multiply(bar_rho_sol, g2)
    
    finite_diff_rhof_x1_sol = np.multiply(np.matmul(rho_f1_sol, div_1)[:-1, :], numerator)
    finite_diff_rhof_x2_sol = np.multiply(np.matmul(div_2,rho_f2_sol)[:, :-1], numerator)
    finite_diff_barrho_x2_sol = np.multiply(np.matmul(div_2,barrho_g2_sol)[:, :-1], numerator)
    div_rho_fg = finite_diff_rhof_x1_sol + finite_diff_rhof_x2_sol + finite_diff_barrho_x2_sol
    
    fig4, axs4 = plt.subplots(1, 1, constrained_layout=True, subplot_kw={"projection": "3d"})
    axs4.set_title('solved div(rho_f+rho_g_u)')
    surf_div = axs4.plot_surface(X_grid[:-1,:-1], Y_grid[:-1, :-1], div_rho_fg, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
    ## visualize the vector field of rho(f+gu)
    fig7, axs7 = plt.subplots(1, 1, constrained_layout=True)
    axs7.set_title('quiver rho(f+gu)')
    axs7.quiver(X_grid[::2,::2], Y_grid[::2, ::2], rho_f1_sol[::2, ::2], (rho_f2_sol + barrho_g2_sol)[::2, ::2], units='width')

    ## visualize solved rho
    fig5, axs5 = plt.subplots(1, 1, constrained_layout=True, subplot_kw={"projection": "3d"})
    surf_div = axs5.plot_surface(X_grid, Y_grid, rho_sol, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    axs5.set_title("solved rho")
    
    ## visualize solved rho_qx
    rho_qx_sol = np.multiply(rho_sol, qx)
    fig8, axs8 = plt.subplots(1, 1, constrained_layout=True, subplot_kw={"projection": "3d"})
    surf_div = axs8.plot_surface(X_grid, Y_grid, rho_qx_sol, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    axs8.set_title("solved rho qx")

    ## visualize solved bar_rho
    fig6, axs6 = plt.subplots(1, 1, constrained_layout=True, subplot_kw={"projection": "3d"})
    surf_div = axs6.plot_surface(X_grid, Y_grid, bar_rho_sol, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    axs6.set_title("solved bar_rho")

    fig3, axs3 = plt.subplots(1, 1, constrained_layout=True)
    axs3.set_title('solved system')
    axs3.quiver(X_grid[::2,::2], Y_grid[::2, ::2], f1[::2, ::2], f2_sol[::2, ::2], units='width')
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


