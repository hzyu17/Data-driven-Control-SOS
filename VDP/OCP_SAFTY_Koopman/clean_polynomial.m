function input_polynomial = clean_polynomial(input_polynomial, tol)
% Wipe out those very small coefficients in a polynomial
issym = false;

if isa(input_polynomial, 'polynomial')
    disp('clean polynomial coefficients')
elseif isa(input_polynomial, 'sym')
    input_polynomial = s2p(input_polynomial);
    issym = true;
    disp('clean polynomial coefficients')
else
    keyboard
end

leng = length(input_polynomial.coefficient(:));
left_coeffs = input_polynomial.coefficient(:);
for i_coef = 1:length(left_coeffs)
    if (abs(left_coeffs(i_coef)) < tol)
        left_coeffs(i_coef) = 0.0;
    end
end
input_polynomial.coefficient = left_coeffs;

if issym
    input_polynomial = p2s(input_polynomial);
end
end

