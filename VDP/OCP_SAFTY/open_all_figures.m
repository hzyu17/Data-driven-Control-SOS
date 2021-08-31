function out = open_all_figures(filename)
% Open all the figure files using one name containing the date in the
% filename, e.g.: 07_03_10_00_Rational_safe_rho_lamb_50000.fig

if contains(filename, 'res')
    date = filename(5:15);
else
    date = filename(1:11);
end
files = {'_Rational_trj_lamb_50000.fig', '_Rational_safe_rho_lamb_50000.fig', '_Rational_safe_divrho_lamb_50000'};
axes = [];
for i_f = 1:length(files)
    i_file = files{i_f};
    full_name = ['res/', date, i_file];
    h_i = openfig(full_name, 'invisible');
    ax_i = gca;
    axes = [axes, ax_i];
end

% Layout of the subfigures
% [    A   |     B  ]
% [    A   |     C  ]

figure('Position', [100, 100, 2048, 1200]);
s1 = subplot(2,2,[1 3]);
grid minor
s2 = subplot(2,2,2);
grid minor
zlim([-10,10]);
s3 = subplot(2,2,4);
grid minor
zlim([-10,10]);
ax1 = axes(1);
ax2 = axes(2);
ax3 = axes(3);
fig1 = get(ax1, 'children');
fig2 = get(ax2, 'children');
fig3 = get(ax3, 'children');

copyobj(fig1, s1);
copyobj(fig2, s2);
copyobj(fig3, s3);

end
