function value = eval_c_over_a(cx, ax, pt)
pvar x1 x2
x = [x1; x2];
dim_cx = length(cx);
value = zeros(dim_cx, 1);
for i_c = 1:dim_cx
    value(i_c) = subs(cx(i_c), x, pt) / subs(ax, x, pt);
end
end

