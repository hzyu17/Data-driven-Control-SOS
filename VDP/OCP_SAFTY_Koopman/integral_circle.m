function integral = integral_circle(input_func, center, radius, x_var, interval)
% To calculate the integration between 2 functions, upper and lower
if isa(input_func, 'polynomial')
    input_func = p2s(input_func);
end
if isa(input_func, 'sym')
    input_func = matlabFunction(input_func, x_var(1));
end

% the lower bound and the upper bound functions
lower_bd = center(2) - sqrt(radius^2 - (x_var(1) - center(1)).^2);
upper_bd = center(2) + sqrt(radius^2 - (x_var(1) - center(1)).^2);

if isa(lower_bd, 'sym')
    lower_bd = matlabFunction(lower_bd);
end
if isa(upper_bd, 'sym')
    upper_bd = matlabFunction(upper_bd);
end

integral = integral2(input_func, interval(1), interval(2), interval(1), interval(2))
% integral = integral2(input_func, interval(1), interval(2), lower_bd, upper_bd)

if integral < 1e-10
    integral = 0.0;
end
end

