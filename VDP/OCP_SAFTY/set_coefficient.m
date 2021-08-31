function input_polynomial = set_coefficient(input_polynomial, coeffs)

poly = 1;
if isa(input_polynomial, 'polynomial')
    disp('clean polynomial coefficients')
elseif isa(input_polynomial, 'sym')
    poly = 0;
    input_polynomial = s2p(input_polynomial);
    disp('clean polynomial coefficients')
else
    keyboard
end

input_polynomial.coefficient(:) = coeffs;
if poly == 0
    input_polynomial = p2s(input_polynomial);
end
end

