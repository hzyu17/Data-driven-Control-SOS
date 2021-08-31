function polynomial = clean_polynomial(polynomial, tol)
% Wipe out small coefficients in a polynomial
% poly_coeffs = poly2basis(polynomial, basis);
% c_div_full = full(poly_coeffs);
% for divi = 1: length(c_div_full)
%     if abs(c_div_full(divi)) < 1e-8
%         c_div_full(divi) = 0.0;
%     end
% end
% poly_coeffs = sparse(c_div_full);
% polynomial  = poly_coeffs' * basis;

leng = length(polynomial.coefficient(:));
left_coeffs = polynomial.coefficient(:);
for i_coef = 1:length(left_coeffs)
    if (abs(left_coeffs(i_coef)) < tol)
        left_coeffs(i_coef) = 0.0;
    end
end
polynomial.coefficient = left_coeffs;

end

