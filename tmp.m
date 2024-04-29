clc; clear;

%%
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));
load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

num_ctl = size(find(subject_index == 1), 1);
num_dep = size(find(subject_index == 2), 1);
num_sui = size(find(subject_index == 3), 1);

% %%
% ctl_id_list = [];
% dep_id_list = [];
% sui_id_list = [];
% 
% for i = 1: length(subject_index)
%     if subject_index(i) == 1
%         ctl_id_list = [ctl_id_list i];
%     elseif subject_index(i) == 2
%         dep_id_list = [dep_id_list i];
%     elseif subject_index(i) == 3
%         sui_id_list = [sui_id_list i];
%     end     
% end
% 
% save('/hd2/research/EEG/data/subject_type_list.mat', ...
%     'ctl_id_list', 'dep_id_list', 'sui_id_list');
    


% suisidal
% sub15 = cur_data(:, :, 15);
% sub23 = cur_data(:, :, 23);
% 
% norm(sub15 - sub23) / norm(sub15)
% 
% %% control
% sub2 = cur_data(:, :, 2);
% 
% norm(sub15 - sub2) / norm(sub15)

%%
clc;

U2 = result(2).U;
sub_coef = U2{3, 1};
labels = [zeros(n_ctl_train, 1); ones(n_sui_train, 1)];
gscatter(sub_coef(:, 1), sub_coef(:, 2), labels);
% scatter(1: 30, sub_coef(:, 1));


