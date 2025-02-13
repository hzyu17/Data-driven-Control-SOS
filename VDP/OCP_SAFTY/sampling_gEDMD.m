function gEDMD_filename = sampling_gEDMD(option)

%% Notes and logs
% 01/12/2021 cellfun, matlabFunction is much more faster than subs() for
% symbolic substitutions
% 08/30/2021 add the data driven part into the whole picture

close all 
clear all

%%
load('polynomials_def')
data_driven_option = option.data_driven_option;
dynamics_option = option.dynamics_option;
% 3. Coefficients of x
c_x = poly2basis(var_poly.x', Psi);

%%
% ============================
% sampling gEDMD for X's and X_dot's
% ============================

u = zeros(dim_m, 1); % input
t = 0; % created but not used 
L = cell(1,dim_m+1);
X = cell(1, dim_m+1);
Y = cell(1, dim_m+1);
for i = 1 : dim_m+1
    Xinit = -5.0 + 10.0 * rand(nx, 10000); % random initial points
    u = zeros(dim_m, 1); 
    if i > 1
        u(i-1, 1) = 1;
    end
    % collecting for the u_i simulation
    Xinit_i = [];
    Yinit_i = [];
    
    %% gEDMD
    % Generate basis Psi and its time derivative expression
%     Psi_sym = p2s(Psi);
    grad_Psi = jacobian(Psi_sym);
    model_dyn = dynamics_definition('sym', dynamics_option);
    model_dyn = model_dyn.F + model_dyn.G*u;
    Psi_dt = grad_Psi * model_dyn;    
    
    % Get initial data Psi(x0)
    data_cell = mat2cell(Xinit, size(Xinit, 1), ones(1, size(Xinit, 2)));
    
    MFun_Psi = matlabFunction(Psi_sym, 'Vars', {var_sym.x});
    Xi = cell2mat(cellfun(MFun_Psi, data_cell, 'UniformOutput', 0));
    
    % Get derivative data on the basis dPsi(x0)
    
    % Important! Use matlabFunction and cellfun to evaluate large
    % scales of data is much more efficient than subs!
    MFun_Psi_dt = matlabFunction(Psi_dt, 'Vars', {var_sym.x});
    dXi =  cell2mat(cellfun(MFun_Psi_dt, data_cell, 'UniformOutput', 0));
    
    % approximation
    M_i = size(Xi, 2);
    A = 1/M_i * dXi * Xi';
    B = 1/M_i * Xi * Xi';
    Li = A * pinv(B);
    figure;
    spy(Li)
    title('Sparsity visualization of the approximated infinitesimal matrix L')
    
    X{1, i} = Xi;
    Y{1,i} = dXi;
    L{1,i} = Li;

end

%% save file
formatOut = 'mm_dd_yy_hh_MM_SS';
ch_date = datestr(now,formatOut);
gEDMD_filename = ['sampling_data/gEDMD_res_', ch_date,'.mat'];
save(gEDMD_filename)