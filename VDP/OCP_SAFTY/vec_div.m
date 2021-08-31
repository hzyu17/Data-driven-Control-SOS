function out = vec_div(func,x_var)
% calculate the divergence of a function vector w.r.t. variable vector
% return sum(d(func(i))/d(x_var(i)))
% Acceptable input function form: f = (f1, ..., fn)'.

% [row, col] = size(func);
% out = [];
% for i_col = 1: col
%     out_irow = 0;
%     for i_row = 1:row
%         out_irow = out_irow + diff(func(i_row, i_col), x_var(i_row));
%     end
%     out = [out, out_irow];
% end
% if col == 1
%     out = out(1);
% end

out = 0;
func_size = size(func);
if length(func) ~= length(x_var)
    disp('not compatible function size and variable size')
    keyboard
end
if min(func_size) > 1
    disp('size not acceptable, only accept vector form functions')
    keyboard
elseif func_size(2) > func_size(1)
    % transpose to size(n,1)
    func = func';
end
% disp('function:')
% vpa(func, 5)
if contains(class(x_var(1)), 'polynomial')
    x_var = p2s(x_var);
end

for i_row = 1:length(func)
    out = out + diff(func(i_row), x_var(i_row));
end
% disp('divergence:')
% vpa(out, 5)
end