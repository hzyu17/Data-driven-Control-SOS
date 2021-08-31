function last_exp_file_name = read_last_data()
% Read from the directory res/ the data of the last experiment
path = [pwd,'/experiments'];
exp_dir = dir(path);
if sum([exp_dir.isdir]) > 0
    exp_dir = exp_dir(~[exp_dir.isdir]);
end
[~, indx] = sort([exp_dir.datenum]);
last_data = exp_dir(indx(end));
last_exp_file_name = [last_data.folder, '/', last_data.name]
end

