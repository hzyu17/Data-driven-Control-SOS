function data_file_name = create_sosprog(option)
%%
% polynomials a, b, c, q, R, w: if it does not exist, please first run polynomials_definition.m
load('polynomials_def.mat')
geometry = domain_definition();
input_regularizer = option.control_penalty_type;
%% sos constraints on OCP equivalent matrix SDP
% add poly_w as decision variables
sos_prog = sosprogram(var_sym.x);

for c_i = 1:dim_m
    c_c_i = c_c(c_i);
    c_c_i = c_c_i{1};
    Qci = Qc(c_i);
    c_c_sym_decl_i = p2s(c_c_i);
    sos_prog = sosdecvar(sos_prog, [c_c_sym_decl_i(1:Qci)]);
end

c_a_sym_decl = p2s(c_a);
sos_prog = sosdecvar(sos_prog, [c_a_sym_decl(1:Qa)]);

%% Base integral
% monomial integrations
lower_bd_X = matlabFunction(geometry.X_cen(2) - sqrt(geometry.r_X^2 - (var_sym.x(1) - geometry.X_cen(1)).^2));
upper_bd_X = matlabFunction(geometry.X_cen(2) + sqrt(geometry.r_X^2 - (var_sym.x(1) - geometry.X_cen(1)).^2));

lower_bd_Xr = matlabFunction(geometry.X_r_cen(2) - sqrt(geometry.r_xr^2 - (var_sym.x(1) - geometry.X_r_cen(1)).^2));
upper_bd_Xr = matlabFunction(geometry.X_r_cen(2) + sqrt(geometry.r_xr^2 - (var_sym.x(1) - geometry.X_r_cen(1)).^2));

integral_vec_a = zeros(Qa, 1);
for li=1:Qa
    cost_fun1_hdl = matlabFunction(Psi_sym(li) .* poly_q_sym ./ (poly_b_sym.^Alph));
%     integral_vec_a(li) = integral_circle(cost_fun1_hdl, geometry.X_cen, geometry.r_X, var_sym.x, [geometry.X_0.xmin, geometry.X_0.xmax]) - ...
%                                               integral_circle(cost_fun1_hdl, geometry.X_r_cen, geometry.r_xr, var_sym.x, [geometry.X_r_coord.xmin, geometry.X_r_coord.xmax]);

% circle
    integral_vec_a(li) = integral2(cost_fun1_hdl, geometry.X_0.xmin, geometry.X_r_coord.xmin, lower_bd_X, upper_bd_X) + ...
                                             integral2(cost_fun1_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, lower_bd_X, lower_bd_Xr) + ...
                                             integral2(cost_fun1_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, upper_bd_Xr, upper_bd_X) + ...
                                             integral2(cost_fun1_hdl, geometry.X_r_coord.xmax, geometry.X_0.xmax, lower_bd_X, upper_bd_X);

% cube
%         integral_vec_a(li) = integral2(cost_fun1_hdl, geometry.X_0.xmin, geometry.X_r_coord.xmin, geometry.X_0.ymin, geometry.X_0.ymax) + ...
%                                              integral2(cost_fun1_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, geometry.X_0.ymin, geometry.X_r_coord.ymin) + ...
%                                              integral2(cost_fun1_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, geometry.X_r_coord.ymax, geometry.X_0.ymax) + ...
%                                              integral2(cost_fun1_hdl, geometry.X_r_coord.xmax, geometry.X_0.xmax,  geometry.X_0.ymin, geometry.X_0.ymax);
end
obj_sos = transpose(c_a_sym(1:Qa)) * integral_vec_a;

if strcmp(input_regularizer, 'L2')
   %% L2 penalty and constraints and objectives on control :w(x), matrix M
    c_w_sym_decl = p2s(c_w);
    sos_prog = sosdecvar(sos_prog, [c_w_sym_decl(1:Qw)]);

    % PSD Matrix M of polynomials
    disp('creating psd matrix M..');
    R_ocp = eye(dim_m);
    size_M = 1+ dim_m;
    M = sym('M', [size_M,size_M]);
    M(1,1) = poly_w_sym;
    for i_c = 2:dim_m+1
        M(1,i_c) = poly_c_sym(i_c-1);
        M(i_c,1) = poly_c_sym(i_c-1);
    end
    M(2:end, 2:end) = R_ocp^-1 .* poly_a_sym; % Should be inv(R) in general cases.
    sos_prog = sosmatrixineq(sos_prog, M, 'Mineq');

    integral_vec_w = zeros(Qw, 1);
    for li=1:Qw
            cost_fun2_hdl = matlabFunction(Psi_sym(li) ./ (poly_b_sym.^Alph));
%             integral_vec_w(li) = integral_circle(cost_fun2_hdl, geometry.X_cen, geometry.r_X, var_sym.x, [geometry.X_0.xmin, geometry.X_0.xmax]) - ...
%                                                       integral_circle(cost_fun2_hdl, geometry.X_r_cen, geometry.r_xr, var_sym.x, [geometry.X_r_coord.xmin, geometry.X_r_coord.xmax]);      

    % circle
%             integral_vec_w(li) = integral2(cost_fun2_hdl, geometry.X_0.xmin, geometry.X_r_coord.xmin, lower_bd_X, upper_bd_X) + ...
%                                                      integral2(cost_fun2_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, lower_bd_X, lower_bd_Xr) + ...
%                                                      integral2(cost_fun2_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, upper_bd_Xr, upper_bd_X) + ...
%                                                      integral2(cost_fun2_hdl, geometry.X_r_coord.xmax, geometry.X_0.xmax, lower_bd_X, upper_bd_X);

    % cube excluding Xr
%             integral_vec_w(li) = integral2(cost_fun2_hdl, geometry.X_0.xmin, geometry.X_r_coord.xmin, geometry.X_0.ymin, geometry.X_0.ymax) + ...
%                                                  integral2(cost_fun2_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, geometry.X_0.ymin, geometry.X_r_coord.ymin) + ...
%                                                  integral2(cost_fun2_hdl, geometry.X_r_coord.xmin, geometry.X_r_coord.xmax, geometry.X_r_coord.ymax, geometry.X_0.ymax) + ...
%                                                  integral2(cost_fun2_hdl, geometry.X_r_coord.xmax, geometry.X_0.xmax,  geometry.X_0.ymin, geometry.X_0.ymax);
    % cube excluding a small set X_excld                                         
    integral_vec_w(li) = integral2(cost_fun2_hdl, geometry.X_0.xmin, geometry.X_excld.xmin, geometry.X_0.ymin, geometry.X_0.ymax) + ...
                                                 integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, geometry.X_0.ymin, geometry.X_excld.ymin) + ...
                                                 integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, geometry.X_excld.ymax, geometry.X_0.ymax) + ...
                                                 integral2(cost_fun2_hdl, geometry.X_excld.xmax, geometry.X_0.xmax,  geometry.X_0.ymin, geometry.X_0.ymax);
    end
    obj_sos = obj_sos + transpose(c_w_sym(1:Qw)) * integral_vec_w;
elseif strcmp(input_regularizer, 'L1')
   %% l1-norm penalty on c(x)
    % vdp and 2d linear
    % beta = 1e3;

    % zero_f_eye_u, vdp
    beta = 5e-4;
    
    for i_indx = 1:dim_m
        Q_cs_i = Qc(i_indx);
        var_name = ['c_l1s', num2str(i_indx)];
        c_l1s_i = mpvar(var_name, [Q_cs_i, 1]);
        c_l1s_i = [c_l1s_i; zeros(nPsi - Q_cs_i, 1)];
        c_l1s_i_sym = p2s(c_l1s_i);
        poly_cl1s_sym_i = transpose(c_l1s_i_sym) * Psi_sym;
        sos_prog = sosdecvar(sos_prog, [c_l1s_i_sym(1:Q_cs_i)]);
        % s(x) >= 0
        sos_prog = sosineq(sos_prog, poly_cl1s_sym_i);
        % s(x) - c(x) >= 0
        sos_prog = sosineq(sos_prog, poly_cl1s_sym_i - poly_c_sym(i_indx));
        % s(x) + c(x) >= 0
        sos_prog = sosineq(sos_prog, poly_cl1s_sym_i + poly_c_sym(i_indx));
        
        % integration and l1-norm objective
        integral_vec_cs = zeros(Q_cs_i, 1);
        for li=1:Q_cs_i
            cost_fun2_hdl = matlabFunction(Psi_sym(li) ./ (poly_b_sym.^Alph));
        % circle geometry.X_excld: a cube
            integral_vec_cs(li) = integral2(cost_fun2_hdl, geometry.X_0.xmin, geometry.X_excld.xmin, lower_bd_X, upper_bd_X) + ...
                                                     integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, lower_bd_X, geometry.X_excld.ymin) + ...
                                                     integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, geometry.X_excld.ymax, upper_bd_X) + ...
                                                     integral2(cost_fun2_hdl, geometry.X_excld.xmax, geometry.X_0.xmax, lower_bd_X, upper_bd_X);
            integral_vec_cs(li) = beta * integral_vec_cs(li);
    
        % cube
        %         integral_vec_cs(li) = integral2(cost_fun2_hdl, geometry.X_0.xmin, geometry.X_excld.xmin, geometry.X_0.ymin, geometry.X_0.ymax) + ...
        %                                              integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, geometry.X_0.ymin, geometry.X_excld.ymin) + ...
        %                                              integral2(cost_fun2_hdl, geometry.X_excld.xmin, geometry.X_excld.xmax, geometry.X_excld.ymax, geometry.X_0.ymax) + ...
        %                                              integral2(cost_fun2_hdl, geometry.X_excld.xmax, geometry.X_0.xmax,  geometry.X_0.ymin, geometry.X_0.ymax);
        end
        obj_sos = obj_sos + transpose(c_l1s_i_sym(1:Q_cs_i)) * integral_vec_cs;
    end
end

%% constraints on the magnitude of u(x): -bdd*a(x) <= c(x) <= bdd*a(x)
% bdd = 1;
% for i_c = 1:dim_m
%     sos_prog = sosineq(sos_prog, bdd * poly_a_sym - poly_c_sym(i_c));
%     sos_prog = sosineq(sos_prog, bdd * poly_a_sym + poly_c_sym(i_c));
% end

 %% saving data
disp('----- Saving solved optimal control data: -----')
date = sprintf('%s', datestr(now,'mm_dd_HH_MM'));
data_driven_type = option.data_driven_option.type;
regularizing_type = option.control_penalty_type;

data_file_name = ['sos_program_data/',date,'_optimal_control_sos_prog_',data_driven_type,'_',regularizing_type,'.mat'];
disp(data_file_name)
save(data_file_name)
end

