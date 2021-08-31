function input_polynomial = clean_polynomial_upper(input_polynomial, tol)
% Wipe out those very large coefficients in a polynomial

if isa(input_polynomial, 'polynomial')
    disp('clean polynomial coefficients')
elseif isa(input_polynomial, 'sym')
    input_polynomial = s2p(input_polynomial);
    disp('clean polynomial coefficients')
else
    keyboard
end

leng = length(input_polynomial.coefficient(:));
left_coeffs = input_polynomial.coefficient(:);
for i_coef = 1:length(left_coeffs)
    if (abs(left_coeffs(i_coef)) > tol)
        left_coeffs(i_coef) = 0.0;
    end
end
input_polynomial.coefficient = left_coeffs;

end

