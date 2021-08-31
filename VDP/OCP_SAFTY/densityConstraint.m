function constraint_sym = densityConstraint(poly_a_sym, poly_b_sym, poly_c_sym, Alph, F, G, nx)
% construct density constraint in sos
Div_F = 0;
Div_G = 0;

syms x1 x2
x = [x1; x2];
for i1 = 1 : nx
    Div_Fi = diff(F(i1), x(i1));
    Div_F = Div_F + Div_Fi;
    Div_Gi = diff(G(i1), x(i1));
    Div_G = Div_G + Div_Gi;
end

% Div_F.coefficient(find(abs(Div_F.coefficient) <= 1e-9)) = 0;
% Div_G.coefficient(find(abs(Div_G.coefficient) <= 1e-9)) = 0;

% F = p2s(F);
% G = p2s(G);
% Div_F = p2s(Div_F);
% Div_G = p2s(Div_G);

poly_ab_sym = poly_a_sym * poly_b_sym;
poly_bc_sym = poly_b_sym * poly_c_sym;

syms x1 x2
grad_a = [diff(poly_a_sym, x1); diff(poly_a_sym, x2)];
grad_ab = [diff(poly_ab_sym, x1); diff(poly_ab_sym, x2)];
grad_bc = [diff(poly_bc_sym, x1); diff(poly_bc_sym, x2)];
grad_c = [diff(poly_c_sym, x1); diff(poly_c_sym, x2)];

term(1) = poly_a_sym * Div_F + grad_a' * F;
term(2) = poly_c_sym * Div_G + grad_c' * G;
term(3) = poly_ab_sym * Div_F + grad_ab' * F;
term(4) = poly_bc_sym * Div_G + grad_bc'* G;

constraint_sym = (1+Alph)*poly_b_sym*sum(term(1:2)) - Alph*sum(term(3:4));
end

