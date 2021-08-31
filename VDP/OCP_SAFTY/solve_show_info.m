function [sos_prog, info] = solve_show_info(sos_prog, option)
%Solve the sos program and print out the solution information

[sos_prog, info] = sossolve(sos_prog, option);

print_sol_info(info)

% disp('------ info.err ------')
% info.err
end

