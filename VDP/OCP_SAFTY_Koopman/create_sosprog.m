function [out] = create_sosprog(poly)
%%
% polynomials a, b, c, q, R, w: if it does not exist, please first run polynomials_definition.m
% close all 
% clear all
% clc

load('polynomials_def.mat')

% domains
%%
geo = domain_definition();
% X = geo.X;
% X_d = geo.X_d;
% X_init = geo.X_init;
% X_r = geo.X_r;
% X_u = geo.X_u;

%%
% dynamics
[F, G] = dynamics_vdp();

if strcmp(class(F), 'polynomial')
    F = p2s(F);
    F = vpa(F, 9);
end

syms x1 x2
x = [x1; x2];

%% sos constraints on OCP equivalent matrix SDP
% add poly_w as decision variables
syms x1 x2
sos_prog = sosprogram([x1, x2]);

for c_i = 1:dim_m
    c_c_sym_i = c_c_sym(c_i);
    c_c_sym_i = c_c_sym_i{1};
    Qci = Qc(c_i);
    sos_prog = sosdecvar(sos_prog, [c_c_sym_i(1:Qci)]);
end

c_a_sym = p2s(c_a);
% poly_a_sym = p2s(poly_a);
sos_prog = sosdecvar(sos_prog, [c_a_sym(1:Qa)]);

% sos constraints on density and a(x)
% sos_prog = sosineq(sos_prog, constraint_sym);
% sos_prog = sosineq(sos_prog, poly_a_sym);

c_w_sym = p2s(c_w);
sos_prog = sosdecvar(sos_prog, [c_w_sym(1:Qw)]);

% PSD Matrix M of polynomials
disp('creating psd matrix M..');
R_ocp = eye(dim_m);
size_M = 1+ dim_m;
M = sym('M', [size_M,size_M]);
M(1,1) = poly_w_sym;
for i_c = 2:dim_m+1
    M(1,i_c) = poly_c_sym(i_c-1);
    M(i_c,1) = poly_c_sym(i_c-1);
end
M(2:end, 2:end) = R_ocp^-1 .* poly_a_sym; % Should be inv(R) in general cases.
M
sos_prog = sosmatrixineq(sos_prog, M, 'quadraticMineq');

%% Base integral
Psi_sym = p2s(Psi);
integral_vec_a = zeros(Qa, 1);

X_excld = geo.X_excld;
X_0 = geo.X_0;

for li=1:Qa
    if poly == 1
        cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym);
        integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_0(4));
    else % rational
        cost_fun1_hdl = matlabFunction(Psi_sym(li) * poly_q_sym / (poly_b_sym^Alph));
        integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + ...
            integral2(cost_fun1_hdl, X_excld(2), X_0(2), X_0(3), X_0(4)) + ...
            integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + ...
            integral2(cost_fun1_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
%         integral_vec_a(li) = integral2(cost_fun1_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun1_hdl, X_excld(2), X_0(2), X_0(3), X_0(4));
%                              + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun1_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    end
end
obj_sos = c_a_sym(1:Qa)' * integral_vec_a;

integral_vec_w = zeros(Qw, 1);
integral_vec_w(1) = (X_0(2) - X_0(1)) * (X_0(4) - X_0(3));
for li=2:Qw
    if poly == 1
        cost_fun2_hdl = matlabFunction(Psi_sym(li), 'Vars', [x1, x2]);
        integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_0(4));
    else % rational
        cost_fun2_hdl = matlabFunction(Psi_sym(li) / (poly_b_sym^Alph));
        integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + ...
            integral2(cost_fun2_hdl, X_excld(2), X_0(2), X_0(3), X_0(4)) + ...
            integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + ...
            integral2(cost_fun2_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
%         integral_vec_w(li) = integral2(cost_fun2_hdl, X_0(1), X_excld(1), X_0(3), X_0(4)) + integral2(cost_fun2_hdl, X_excld(2), X_0(2), X_0(3), X_0(4));
%                              + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_0(3), X_excld(3)) + integral2(cost_fun2_hdl, X_0(1), X_0(2), X_excld(4), X_0(4));
    end
end
obj_sos = obj_sos + c_w_sym(1:Qw)' * integral_vec_w;


%% sos program for safety verification
load('safety_prog.mat')

%% save data
save('sos_prog')
end

