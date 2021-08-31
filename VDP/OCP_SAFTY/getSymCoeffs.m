function c_out = getSymCoeffs(poly_sym, x, monomial_basis)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% poly_flag = 0;
% if isa(poly_sym, 'polynomial')
%     poly_sym = p2s(poly_sym);
%     poly_flag = 1;
% end
[coeff_c_sym, T] = coeffs(poly_sym, x);
coeff_c_sym = [coeff_c_sym];
total_coeffs = length(monomial_basis);
c_out = zeros(total_coeffs,1);
for i = 1:length(T)
    monomial_i = T(i);
    coeff_monomi = s2p(monomial_i);
    term_i = poly2basis(coeff_monomi, monomial_basis) .* s2p(coeff_c_sym(i));
    c_out = c_out+term_i;
end
% if poly_flag == 1
%     c_out = s2p(c_out);
% end
end

