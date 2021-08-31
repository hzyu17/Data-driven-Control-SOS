function out = vec_grad(func,x_var)
out = [];
for ii = 1:length(x_var)
    out = [out; diff(func, x_var(ii))];
end

end

