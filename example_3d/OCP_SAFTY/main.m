clear all
close all
clc

%%
poly_flg = 0; %polynomial or rational

coeff = 5e12;
polynomials_definition(coeff);
% create_safety_program();
create_sosprog(poly_flg);

%% experiments 
% lambdas = [1e40, 1e50, 1e60, 1e70, 1e72, 1e73, 1e74, 1e75, 1e76, 1e77, 1e80];
% for higher order rationals:
% lambdas = [1e110, 1e112, 1e114, 1e116, 1e118,1e120,1e122,1e124,1e126,1e128, 1e129];
lambdas = [1e6 1e8];
epss = [0.1, 0.2, 0.3, 0.4];% for polynomials:

geo = domain_definition();
eps = 0;

% for lamda = lambdas(1:end)
for lamda = 1e5
%     for eps = epss(1:end)
%     for eps = 0.2
    disp('============ Lambda ============')
    lamda
    
    if poly_flg == 1
        ocp_safety_poly_main(lamda, eps, geo);
    else
        data_file_name = ocp_safety_rational_main(lamda, eps, geo);
    end
    %% Drawings
    plot_result(lamda, eps, poly_flg, geo, data_file_name)
    aa = 1;
%     end
end

%% see figs
% lbd = 0;
% ee = 0;
% if poly_flg == 1
%     rat_poly = 'Poly';
% else
%     rat_poly = 'Rational';
% end
% 
% name_div = ['res/',rat_poly,'_divrho_lamb_',num2str(lbd),'_eps_',num2str(ee),'.fig'];
% name_rho = ['res/',rat_poly,'_rho_lamb_',num2str(lbd),'_eps_',num2str(ee),'.fig'];
% name_trj_fig = ['res/',rat_poly,'_trj_lamb_',num2str(lbd),'_eps_',num2str(eps),'.fig'];
% name_safety_rho = ['res/',rat_poly,'_safe_rho_','lamb_',num2str(lbd),'_eps_',num2str(eps),'.fig'];
% name_safety_div = ['res/Rational_safe_divrho_','lamb_',num2str(lbd),'_eps_',num2str(eps),'.fig'];
% 
% openfig(name_div,'visible');
% openfig(name_rho,'visible');
% openfig(name_trj_fig,'visible');
% openfig(name_safety_rho,'visible');
% openfig(name_safety_div,'visible');

%% comparison figures
% fig_handles = [];
% fig_axis = [];
% % lambdas = [1 2 3 4 5 6 7 8 9];
% lambdas = [.1 .2 .3 .4 .5 .6 .7 .8 .9];
% ee = 0;
% if poly_flg == 1
%     rat_poly = 'Poly';
% else
%     rat_poly = 'Rational';
% end
% for lbd = lambdas(1:end)
%     name_trj = ['res/',rat_poly,'_trj_','lamb_',num2str(lbd),'_eps_',num2str(ee),'.fig'];
%     h_i = openfig(name_trj, 'reuse');
%     axis_i = gca;
%     fig_axis = [fig_axis; axis_i];
%     fig_handles = [fig_handles; h_i];
% end
% 
% h_together = figure;
% subplots_t = tight_subplot(3, 3, [0.01 0.01]);
% for ii = 1:3
%     for jj = 1:3
%         indx = (ii-1)*3+jj;
%         s_ij = subplots_t(indx);
%         fig_ij = get(fig_axis(indx), 'children');
%         copyobj(fig_ij, s_ij);
%     end
% end