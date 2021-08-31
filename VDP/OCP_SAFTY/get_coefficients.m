function [out] = get_coefficients(input_polynomial)
%% get the coefficient of a polynomial using a symbolic input 

if isa(input_polynomial, 'polynomial')
    disp('clean polynomial coefficients')
elseif isa(input_polynomial, 'sym')
    input_polynomial = s2p(input_polynomial);
    disp('clean polynomial coefficients')
else
    keyboard
end

out = input_polynomial.coefficient(:);

end


