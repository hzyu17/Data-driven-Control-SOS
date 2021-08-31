function output = print_sol_info(info)
% print the solution info given by the sdp solver

disp('==================Solution check for the optimal control problem================')
disp('----------- feasibility check --------------')
disp('info.pinf')
info.pinf
disp('info.dinf')
info.dinf
if info.pinf == 1
    disp('Primal infeasibility')
elseif info.dinf == 1
    disp('Dual infeasibility')
else
    disp('------------- Found feasible optimal control rational function ---------------')
end

disp('------- numerical issue -----------')
info.numerr
if info.numerr == 2
    disp('Attention! numerical issues! Do not trust the SDP results')
    keyboard;
elseif info.numerr == 1
    disp('Attention! numerical inaccuracy!')
else
    disp('NO numerical issues, great!')
end

disp('------- info.feasratio -------')
info.feasratio

output = 0;

end

