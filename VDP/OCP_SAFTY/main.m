%% Optimal control for nonlinear dynamical systems with L2 cost and safty constraints using SOS
% assume we know the dynamics
% process: 
% OCP Part:
% 1. Read the dynamics f(x) from functions or Koopman
% 2. Read polynomials and Create SOS program
% 2. Define the reward matrix q(x), R(x), and calculate integral OCP cost
% 3. Define the Matrix equality and SOS feasibility problem
% 4. Add the density SOS constraints

% Safty Part:
% 1. Read the initial sets, unsafe sets, and reachable sets 
% 2. Formulate the SOS constraints w.r.t. X_0, X_u and X_r
% 3. Solving the whole program

% Log
% 02/20/2021: Start from OCP problem: recover the polynomials without
%                             Koopman
% 02/21/2021: Formulate and solve the OCP problem in matrix SDP in SOSTOOLS.
% 03/07/2021: Formulate the safety cost in the formulation
% 06/17/2021: fixed the integration area bug and changed the Matrix SOS
%                             criteria from '’quadraticMineq’' to 'Mineq’, added the multiple unsafe
%                             set case, and trimed the objective function
%                             decimal before putting into the sos solver
% 06/20/2021: found out in the paper what is 'dual infeasible', and test if
%                             it is brought up by the negative in the objective function: it's not
%                             brought up by the negativity.
% 06/21/2021: do we need to convert the lambda to symbolic? tested with
%                             this possibility. Also changed the constant coefficient of a(x) to 1, try
%                             with it: Not working. Still the dual
%                             infeasibility problem. The dual infeasibility
%                             criteria is in the line 690 in the sedumi.m
%                             check that.
% 06/22/2021: The objective function has some very big terms (1e3) and
%                            also some very small terms (1e-13). The small
%                            terms may comes from the two integration
%                            difference of integrate over (X_0 - X_excld)
%                            and (X_r - X_excld). Try directly with
%                            (X_r-X_excld): helped a little, the rho_safe
%                            seems better. (when lambda==0)
% 06/22/2021: test truncation of the objective function before sending into
%                             the solver: nearly the same with the previous
%                             one. 
% 06/22/2021: need to check the stop criteria of the algorithm in the
%                             sedumi.m.
% 06/22/2021: tested the real number assumption on the 
% 06/30/2021: removed the truncations in the objective setting part, and it
%                             gives the feasible solution in the case of feasibility problem
% 07/01/2021: the integration is symmetric, so all the polynomials with odd
%                             degree will be zero. will that pose problems? try with isymetric settings
% 07/04/2021: rescaled the objective function by dividing the min
%                             coefficient, to provide a better numerical condition; changed the
%                             integration from the box to circle.
%07/04/2021: added a function to open all the related figures given one of
%                             names. including the trajectory, the rho_safe, and the div(rho_safe(f+gu))
%07/12/2021: re-organized the structure of the code: separate the
%                            pre-plotting and the plotting action, separated the optimal control and
%                            the safety verification problem.

clear all
close all
clc

%%
% poly_flg = 0; %polynomial or rational
dynamics_option = 'vdp';
data_driven_option.type = 'data_driven'; % choice between 'data_driven' and 'model_based'
data_driven_option.approx = 'gEDMD';
%% polynomial definitions
polynomials_definition(dynamics_option); % availables: linear, minus_cubic, rantzer_example, zero_f_eye_u, vdp

%% gEDMD
if strcmp(data_driven_option.type, 'data_driven')
    gEDMD_filename = sampling_gEDMD(data_driven_option);
else
    gEDMD_filename = '';
end

%% solve ocp
% for lamda = 1e5 % corresponding to van der pol dynamics
% lambda = 50; % 2-dim integrator
% for lambda = [0, 100, 200, 300, 400, 1000] % L1 vdp
for lambda = [10, 20, 30, 40, 50] % L1 vdp data-driven
% for lambda = [10, 20, 30, 40, 50] % L2 vdp
    sos_prog_file = create_sosprog('L1');
    ocp_data_file_name = solve_optimal_control(lambda, sos_prog_file, data_driven_option, gEDMD_filename);

    % plotting ocp results
    % ocp_data_file_name = 'experiments/07_02_00_03_optimal_control_results_SOS_Rational_lbd_50000.mat';

    plot_ocp_results(ocp_data_file_name);
end

%% solve safety
% ocp_data_file_name = 'experiments/07_02_00_03_optimal_control_results_SOS_Rational_lbd_50000.mat';
% ocp_data_file_name = 'experiments/07_02_16_56_optimal_control_results_SOS_Rational_lbd_50000.mat';
% ocp_data_file_name =
% 'experiments/07_02_23_39_optimal_control_results_SOS_Rational_lbd_50000.mat';

% safety_prog_file = create_safety_program();
%%
safety_prog_file = '';
data_file_name = solve_safety_verification(ocp_data_file_name, safety_prog_file);
% ----------- plotting ocp+safety -------------
% load('experiments/07_05_00_53_results_SOS_Rational_lbd_10000.mat')
plot_from_file(data_file_name)

%% load the last time experiment data and draw
% last_time_experiment_file_name = read_last_data();
% load(last_time_experiment_file_name);
% plot_result(lambda, eps, poly_flg, geo, data_file_name)

%%
% data_file_name = ['experiments/', '07_08_15_14_results_SOS_Rational_lbd_0_ocp_feasratio_-0.69971_safety_feasratio_0.96015.mat'];
% print_info_from_file(data_file_name);
