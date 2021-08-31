function [outputArg1,outputArg2] = change_figure_size(input_figure)
fig = open(input_figure)
xlim([-2.5, 1.5])
ylim([-0.9, 3.1])
end

