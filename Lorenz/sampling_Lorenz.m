%% sampling data points X_i's and Y_i's

nx = 3; % variable dimension
m = 1; % input dimension
ns = 1e4; % number of simulations

%% linearized dynamics, sampling directly
T = 9990;
dt = 1e-2;
% [X, Y] = VanDerPol_one_pt(T, dt);
% [X, Y] = Inverted_Pendulum(T, dt);

%% exact model dynamics, sampling after solving the system ODE
rel_tol = 1e-6;
abs_tol = 1e-6;
absxlim = [5.01, 5.01, 5.01];
tspan = [0: dt: 2];
odeoptions = odeset('RelTol',  rel_tol, 'AbsTol', abs_tol, 'Events', @(t,x) event_out_of_limit(t,x,absxlim));

% generate random initial points and collector for the next point

% - Shape: X and Y are cells of simulation data.

Xinit = -5.0 + 10.0 * rand(nx, 10000); % random initial points
u = 0; % input

for i = 1 : m+1
    if i == 2
        u = 1;
    end
    % collecting for the u_i simulation
    Xi = [];
    Yi = [];
    % solve ODE and sample
    for si = 1 : ns
        xinit = Xinit(:, si);
        [t,x,xe,ye,ie] = ode45(@(t,x) model_dynamics_Lorenz(t,x,u), tspan, xinit, odeoptions);
        
        if isempty(ie) % abs(x) within absxlim
            t = t(1:2);
            x = x(1:2, :);
        else % x exceeds absxlim
            if length(t) >= 3
                t = t(1:2);
                x = x(1:2, :);
            else % the first iteration exceeds absxlim
                continue
            end
        end
        Xi = [Xi, x(1, :)'];
        Yi = [Yi, x(2, :)'];
    end
    
    if i==1
        X = Xi;
        Y = Yi;
    else
        X = {X,Xi};
        Y = {Y,Yi};
    end
    
end

save('sampling_data.mat')