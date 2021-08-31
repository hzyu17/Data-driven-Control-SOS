function title = figure_title(dynamics_option, lambda, ocp_info, sos_safety_info)
% title for the figures
    if isempty(sos_safety_info)
        title = [dynamics_option, '_lambda:', num2str(lambda),'_ocp feasratio:',num2str(ocp_info.feasratio)];
    else
        title = [dynamics_option, '_lambda:', num2str(lambda),'_ocp feasratio:',num2str(ocp_info.feasratio),'_safety feasratio:',num2str(sos_safety_info.feasratio)];
    end
end
