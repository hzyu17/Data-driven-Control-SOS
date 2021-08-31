% % Test on l1 norm objective SOS optimizaitons
% Author: Hongzhe Yu
% Objective: Test the SOS optimization with the l1 norm on its polynomial
% coefficients.
% --------------------------------------------------------------------------------
clc
clear 
pvar x

% Monomial
order = 4;
Psi = monomials(x, 0:order);
nPsi = size(Psi,1);

% Polynomial
c_f = mpvar('c', nPsi,1);
poly_f = p2s(c_f'*Psi);
c_f_sym = p2s(c_f);

% SOS program
x = p2s(x);
prog = sosprogram(x);
prog = sosdecvar(prog, c_f_sym);

% constraints on df
constraint_df = diff(poly_f, x);
prog = sosineq(prog, constraint_df-2);
% soseq() sosineq()
% prog = sosineq(prog, poly_f-2);

% % information print out
% symdecvartable = prog.symdecvartable
% constraint_df - 2
poly_f
At = full(prog.expr.At{1})
b = full(prog.expr.b{1})
Z = full(prog.expr.Z{1})
% constr = constraint_df - 2

% Solve
[prog, info] = sossolve(prog);

% Get result
c_f_sol = sosgetsol(prog, c_f_sym)
poly_f_sol = sosgetsol(prog, poly_f)

%%
% ---------- Program 2: with l1-norm objective --------------

% SOS program
prog2 = sosprogram(x);
prog2 = sosdecvar(prog2, c_f_sym);

% constraints on df
constraint_df = diff(poly_f, x);
prog2 = sosineq(prog2, constraint_df-2);
% prog2 = sosineq(prog2, poly_f);

% %% l1-norm minimization
% % l1-norm variables
% c_abs_sym = p2s(mpvar('d', nPsi,1));
% prog2 = sosdecvar(prog2, c_abs_sym);
% 
% % l1 constraints
% for i_var = 1:size(c_abs_sym,1)
%     prog2 = sosineq(prog2, c_abs_sym(i_var) - c_f_sym(i_var));
%     prog2 = sosineq(prog2, c_abs_sym(i_var) + c_f_sym(i_var));
% end
% 
% % l1 objective
% obj = ones(size(c_abs_sym,1),1)'*c_abs_sym;
% sossetobj(prog2, obj);

% % information print out
% symdecvartable = prog2.symdecvartable
% constraint_df - 2
% At = full(prog2.expr.At{1})
% b = full(prog2.expr.b{1})
% Z = full(prog2.expr.Z{1})
% constr = constraint_df - 2

% Solve
[prog2, info] = sossolve_l1_min(prog2);

% Get result
c_f_sol2 = sosgetsol(prog2, c_f_sym)
% c_f_abs_sol2 = sosgetsol(prog2, c_abs_sym)
poly_f_sol2 = sosgetsol(prog2, poly_f)

%% Verify the two results
coeff_f_sol = getSymCoeffs(poly_f_sol, x, Psi);
l1_f_sol = sum(abs(coeff_f_sol))

coeff_f_sol2 = getSymCoeffs(poly_f_sol2, x, Psi);
l1_f_sol2 = sum(abs(coeff_f_sol2))

%% Drawings of the results
figure
x_pts = linspace(-50,50);
y_pts = subs(poly_f_sol, x_pts);
plot(x_pts, y_pts)
grid minor

figure
x_pts = linspace(-50,50);
y_pts = subs(poly_f_sol2, x_pts);
plot(x_pts, y_pts)
grid minor

