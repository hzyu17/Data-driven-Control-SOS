% Van der Pol Oscillator: d_x1 = x2; d_x2 = (1-x1^2)x2 - x1 + u
clear all
close all
clc

%% collect data points X_i's and Y_i's
T = 9990;
dt = 0.01;
[X, Y] = VanDerPol_one_pt(T, dt);
% [X, Y] = Inverted_Pendulum(T, dt);

%% polynomial generation
clc
n = 2; % variable dimension
m = 1; % input dimension
pvar x1 x2 u
x = [x1;x2];

% monomials
highest_order = 6;
Phix = monomials(x, 0:highest_order);

% polynomials a, b, alpha
poly_b = x'*x;
alpha = 4;
poly_a = 1;

% polynomial c(x) is the one we want to get. 
q_c = 4;
Qc = nchoosek(q_c + n, n);
total_coeffs = nchoosek(highest_order+n, n);
c_c = mpvar('c', [Qc, 1]);
c_c(1) = 0;
z = zeros(total_coeffs-Qc, 1);
p = polynomial(z);
c_c = [c_c;p];

poly_c = c_c' * Phix;

% coefficients
c_b = poly2basis(poly_b, Phix);
c_a = poly2basis(poly_a, Phix);
c_ab = poly2basis(poly_a * poly_b, Phix);
c_bc = c_c .* poly_b;
c_x = poly2basis(x', Phix);

% x = c_x' * Phix;

% the final variable: all the polynomial coefficient combined
d = [c_a; c_c];

%% Koopman approximation
% get A and B matrices
A = {};
B = {};
Xis = {};
Yis = {};

for i = 1:m+1
    Xi = X(2*i-1 : 2*i, :);
    Yi = Y(2*i-1 : 2*i, :);
    Ti = length(Xi);
    
    sum_XX = zeros(Phix.matdim(1), Phix.matdim(1));
    sum_XY = zeros(Phix.matdim(1), Phix.matdim(1));

    for iter = 1:Ti
        x_step = Xi(1:2, iter);
        y_step = Yi(1:2, iter);
        px_step = subs(Phix, x, x_step);
        px_step = px_step.coefficient';
        py_step = subs(Phix, x, y_step);
        py_step = py_step.coefficient';
        sum_XX = sum_XX + px_step * px_step';
        sum_XY = sum_XY + py_step * px_step';
    end
    
    if length(A) <= 0
        A = sum_XX ./ Ti;
        B = sum_XY ./ Ti;
        Xis = Xi;
        Yis = Yi;
    else
        A =  {A, sum_XX ./ Ti};
        B =  {B, sum_XY ./ Ti};
        Xis = {Xis, Xi};
        Yis = {Yis, Yi};
    end
end

%% calculate Ki matrices form an argmin optimization
clc
K = {};
for i = 1:m+1
    A_pinv = pinv(A{i});
    if i == 1
        K = B{i} * A_pinv;
    else
        K = {K, B{i} * A_pinv};
    end
end

% clean up K, make the small numbers to zero
% for i = 1:m+1
%     K_i = K{i};
%     for j=1:length(K_i(:,1))
%         for k=1:length(K_i(1,:))
%             if abs(K_i(j,k)) < 1e-8
%                 K{i}(j,k) = 0.0;
%             end
%         end
%     end
% end

%K are right
% c_x' * K * Phix = [y1; y2] 
% in the Van der Pol dynamical system should be:
% [x1+x2*dt; x2+[(1-x1^2)x2 - x1]*dt]: [x1+0.01x2; 1.01x2-0.01x1-0.01x1^2*x2]
% c_x' * K{1} * Phix
% [x1+0.01x2; 1.01x2-0.01x1-0.01x1^2*x2+0.01]
% c_x' * K{2} * Phix

%% construct L0, L1
clc
% L
for i=1:m+1
    if i == 1
        L = (K{1} - eye(total_coeffs)) ./ dt;
    else
        L = {L, (K{i} - K{1}) ./ dt;};
    end
end

% clean L matrix
% for i =1:2
%     for j = 1:length(L{i}(:,1))
%         for k = 1:length(L{i}(1,:))
%             if L{i}(j,k) < 1e-8
%                 L{i}(j,k) = 0.0;
%             end
%         end
%     end
% end

% [x2; (1-x1^2)x2-x1]
% c_x' * L{1} * Phix
% [0; 1]
% c_x' * L{2} * Phix

pvar x1 x2

% P
for j = 1:m+1
    poly_div = c_x' * L{j} * Phix;

    % j == 1: div(F) = div([K0(x1), ..., K0(xn)]) = div(c_x' * L0 * Phix)
    % j >=2: div(Gj) = div(c_x' * Lj * Phix)
    div = diff(poly_div(1, 1), x1) + diff(poly_div(2, 1), x2);
    
    % polynomial coefficients of this divergence term
    % c_div' * Phix is the div(F) or div(Gj) for [x1, x2]
    
    % j=1: 1-x1^2
    % j=2: 0
    div = clean_polynomial(div)
    
    if j == 1
        P = L{1} + div .* eye(total_coeffs);
    else
        P = {P, L{j} + div .* eye(total_coeffs)};
    end
end

%% formulate SOS optimization
% SOS constraint
% clc
% for j = 0:m
%     if j == 0
%         continue
%     end
%     if j == 1
%         sum_cc = P{j} * Phix;
%         sum_bc = P{j} * Phix;
%     else
%         sum_cc = sum_cc + P{j+1} * Phix;
%         sum_bc = sum_bc + P{j+1} * Phix;
%     end
%     
% end

%%
clc
% constraints
constraint = (1+alpha) .* poly_b .* (c_a' * (P{1} * Phix) + c_c' * (P{2} * Phix)) - alpha .* (c_ab' * (P{1} * Phix) + c_bc' * (P{2} * Phix)) >= 0;

%% Try a simple case with SOSTOOLS to verify the symbolic toolbox
% clc 
% syms x1 x2;
% program1 = sosprogram([x1;x2]);
% c_c_sym = p2s(c_c);
% c_c_sym(2)
% program1 = sosdecvar(program1, [c_c_sym(2:end)]);
% 
% p_temp = c_c' * (Phix.*Phix);
% p_temp_sym = p2s(p_temp);
% 
% program1 = sosineq(program1, p_temp_sym);
% program1 = sossolve(program1);
% sol_p = sosgetsol(program1, p_temp_sym)
% [Q, Z] = findsos(sol_p)
% sol_c_c = sosgetsol(program1, c_c_sym)

%% Use SOSTOOLS to deal with the SOS problem in the data-driven project
clc 
% =========== declare SOS program and independent polynomial variables
syms x1 x2;
program2 = sosprogram([x1;x2]);

% =========== declare decision vairables: the c's
c_c_sym = p2s(c_c);
program2 = sosdecvar(program2, [c_c_sym(2:Qc)]);

poly_c_sym = p2s(poly_c);

% =========== declare constraints: the density SOS constraint

poly_div_gbc = c_bc' * P{2} * Phix;
poly_div_gbc_sym = p2s(poly_div_gbc);

poly_div_gc = c_c' * P{2} * Phix;
poly_div_gc_sym = p2s(poly_div_gc);

poly_div_fa = c_a' * P{1} * Phix;
poly_div_fa_sym = p2s(poly_div_fa);

poly_div_fab = c_ab' * P{1} * Phix;
poly_div_fab_sym = p2s(poly_div_fab);

poly_b_sym = p2s(poly_b);

constraint_sym = (1+alpha) * poly_b_sym * (poly_div_fa_sym + poly_div_gc_sym) - alpha * (poly_div_fab_sym + poly_div_gbc_sym);

% poly_main = (1+alpha) * poly_b * (c_a' * (P{1} * Phix) + c_c' * sum_cc) - alpha * (c_ab' * (P{1} * Phix) + c_bc' * sum_bc);
% poly_main_sym = p2s(poly_main);
% =========== declare objective function 
objective1 = 0;
for i = 2:15
    i_c = c_c_sym(i, 1);
    objective1 = objective1 + i_c;
end

program2 = sossetobj(program2, objective1);

% poly_temp = clean_polynomial(poly_temp)poly_main_sym = p2s(poly_main);
program2 = sosineq(program2, constraint_sym);

% =========== solve the program
program2 = sossolve(program2);

% =========== get SOS solution
sol_poly_main_sym = sosgetsol(program2, constraint_sym);
sol_poly_main_p = s2p(sol_poly_main_sym);

sol_poly_main_p.coefficient;
sol_c_c_sym = sosgetsol(program2, c_c_sym)

%% plot system under the calculted input and inittial state
% dynamics
pvar x1 x2
sol_poly_c_sym = sosgetsol(program2, poly_c_sym);
sol_poly_c = s2p(sol_poly_c_sym);

sol_poly_c = sol_poly_c/200.0;
sol_poly_c = clean_polynomial(sol_poly_c)
u_solved = sol_poly_c;
% u_solved = 0.9015*x1^2*x2 + 0.0251*x2^3 + 0.0241*x2^2 - 1.2505*x2;
x1dot = x2;
x2dot = (1-x1^2)*x2 - x1 + u_solved;
xdot = [x1dot; x2dot];

x_pts = [];
initial_x = [0.01, -0.01, 0.01, -0.01, 0.02, -0.02, 0.02, -0.02, 0.0, 0.0, 0.1, -0.1];
initial_y = [0.01, -0.01, -0.01, 0.01, 0.02, -0.02, -0.02, 0.02, 0.1, -0.1, 0.0, 0.0];

initial_x = [2.5, 2.5, -2.5, -2.5, 0.02, -0.02, 0.02, -0.02];
initial_y = [2.5, -2.5, 2.5, -2.5, 0.02, -0.02, -0.02, 0.02];

figure
hold on

% for kk = 1:4
    for i=1:2000
        if i == 1
%             x_pts = [initial_x(kk); initial_y(kk)];
            x_pts = [-5; -5];
        else
            u_solved_value = subs(u_solved, x, x_pts(1:2, i-1))
            x_dot_poly_value = subs(xdot, x, x_pts(1:2, i-1));
            x_dot_pt = x_dot_poly_value.coefficient';
            x_pts(1:2, i) = x_pts(1:2, i-1) + x_dot_pt.*dt;
            plot(x_pts(1,i),x_pts(2,i),'.')
        end
    end
% end


% example feasibility problem solving
% f1 = c_c' * Phix >= 0;
% [info, dopt, sossol] = sosopt(f1, x, []);




