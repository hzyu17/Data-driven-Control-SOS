function [X,Y] = VanDerPol_one_pt(T, dt)
% Van der Pol Oscillator: d_x1 = x2; d_x2 = (1-x1^2)x2 - x1 + u

%% collect data points X_i's and Y_i's

% Shape: X and Y are in shape: [(m+1)*2, T/dt]. Every 2 rows corresponds to
% an input from 0 to m. the number of columns corresponds to the number of
% data points in one simulation.

n = 2; % variable dimension
m = 1; % input dimension

% input
u = 0;

Xi = -5.0 + 10.0 * rand(2, 10000);
Yi = zeros(2, 10000);    
for i = 1:m+1
    % dynamics
    count = 1;
    
    if i == 2
        u = 1;
    end
    for x_step = Xi(1:2, :)

        vel = [x_step(2,1); (1-x_step(1,1)^2)* x_step(2,1) - x_step(1,1) + u];

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
