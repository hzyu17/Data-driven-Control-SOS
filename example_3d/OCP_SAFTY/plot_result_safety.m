function plot_result_safety(X_0, X_init, X_u, X_r)
% Plotting of the safety constraint problem
step = 0.05;
x1_linspace = X_0(1):step:X_0(2);
x2_linspace = X_0(3):step:X_0(4);

figure
grid minor
[X1, Y1] = meshgrid(x1_linspace, x2_linspace);
[X2, Y2] = meshgrid(x1_linspace, x2_linspace);
[X3, Y3] = meshgrid(x1_linspace, x2_linspace);

% X_init
outside_u = find(feval(matlabFunction(X_init), X1, Y1)<0);
X1(outside_u) = nan; 
Y1(outside_u) = nan;

% X_u
outside_u = find(feval(matlabFunction(X_u), X2, Y2)<0);
X2(outside_u) = nan; 
Y2(outside_u) = nan;

% X_r
outside_r = find(feval(matlabFunction(X_r), X3, Y3)<0);

X3(outside_r) = nan;
Y3(outside_r) = nan;

hold on
plot(X1, Y1, 'Marker', '.', 'Color', 'g')
plot(X2, Y2, 'Marker', '.', 'Color', 'b')
plot(X3, Y3, 'Marker', '.', 'Color', 'r')

xlim(X_0(1:2))
ylim(X_0(3:4))
end

