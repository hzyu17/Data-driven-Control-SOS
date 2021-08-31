function output = print_info_from_file(data_file)
%% load data file
load(data_file)
disp('dynamics: ')
dynamics_option = 'rantzer_example'
save(data_file)

disp('====== optimal control problem ======')
disp('------ polynomial definitions ------')
total_deg
deg_a
deg_c

disp('------ objective function -------')
obj_sos

disp('------ polynomial solution ------')
ax_sol
cx_sol
bx_sol

if exist('ocp_info') == 1
    disp('------ solution info ------')
    print_sol_info(ocp_info)
end


disp('====== verification problem ======')
disp('------ dynamics ------')
if exist('dynamics') == 1
    dynamics.F
    dynamics.G
else
    F
    G
end
disp('------ polynomials ------')
rho_safe_sol

if exist('sos_safety_info') == 1
    disp('------ solution info safety ------')
    print_sol_info(sos_safety_info)
end

%% plotting
plot_from_file(data_file);

end

