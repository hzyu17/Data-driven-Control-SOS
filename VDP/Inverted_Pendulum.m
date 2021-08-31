function [X, Y] = Inverted_Pendulum(T, dt)
%Inverted pendulum dynamical system: 
% x1_dot = x2, x2_dot = sinx1 - 0.5x2 + u

%% collect data points X_i's and Y_i's

% Shape: X and Y are in shape: [(m+1)*2, T/dt]. Every 2 rows corresponds to
% an input from 0 to m. the number of columns corresponds to the number of
% data points in one simulation.

n = 2; % variable dimension
m = 1; % input dimension

% input
u = 0;

% collecting for each dimension in the input
% x1 = linspace(-5.0, 5.0, 100);
% x2 = linspace(-5.0, 5.0, 100);
% 
% count_x = 1;
% Xi = zeros(2, length(x1) * length(x2));
% Yi = zeros(2, length(x1) * length(x2));
% for i_x1 = x1(1, :)
%     for i_x2 = x2(1, :)
%         Xi(:, count_x) = [i_x1; i_x2];
%         count_x = count_x + 1;
%     end
% end

Xi = -pi + 2 * pi * rand(2, 10000);
Yi = zeros(2, 10000);    
    
for i = 1:m+1
    % dynamics
    count = 1;
    
    if i == 2
        u = 1;
    end
    for x_step = Xi(1:2, :)
                
        d_x1 = x_step(2,1);
        d_x2 = sin(x_step(1,1)) - 0.5*x_step(2,1) + u;

        vel = [d_x1; d_x2];

        Yi(1:2, count) = x_step + vel .* dt;
        count = count + 1;
    end
    if i == 1
        X = Xi;
        Y = Yi;
    else
        X = [X; Xi];
        Y = [Y; Yi];
    end
    
end

end

