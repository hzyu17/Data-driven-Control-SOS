function out = vec_div(func,x_var)
out = 0;
for i_x = 1: length(x_var)
    out  = out + diff(func(i_x), x_var(i_x));
end
end

