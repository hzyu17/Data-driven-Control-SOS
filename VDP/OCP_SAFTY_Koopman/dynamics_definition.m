function dynamics = dynamics_definition(var_type, dynamics_type)

if strcmp(var_type, 'sym')
    syms x1 x2
    x = [x1;x2];
elseif strcmp(var_type, 'polynomial')
    pvar x1 x2
    x = [x1;x2];
end

switch dynamics_type
    case 'vdp'
        % van der pol
        F = [x(2); (1 - x(1)^2) * x(2) - x(1)];
        G = [0; 1];
    case 'zero_f_eye_u'
        % F = 0, G=I
        F = [0; 0];
        G = [1, 0; 0, 1];
    case 'linear'
        % linear dynamics
        F = [-x(1) ; - x(2)];
        G = [1, 0; 0, 1];
    case 'minus_cubic'
        % testing dynamics
        F = [-x(1) - x(2); x(1) - x(2)^3];
        G = [1, 0; 0, 1];
    case 'rantzer_example'
        % dynamics in Rantzer's paper
        F = [x(2); -x(1) + x(1)^3./3 - x(2)];
        G = [1, 0; 0, 1];
end

dynamics.F = F;
dynamics.G = G;

end

