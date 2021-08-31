function out = vec_div(func,x_var)
[row, col] = size(func);
out = [];
for i_col = 1: col
    out_irow = 0;
    for i_row = 1:row
        out_irow = out_irow + diff(func(i_row, i_col), x_var(i_row));
    end
    out = [out, out_irow];
end
if col == 1
    out = out(1);
end
end