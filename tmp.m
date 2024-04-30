clc; clear;

%%
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));
load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

num_ctl = size(find(subject_index == 1), 1);
num_dep = size(find(subject_index == 2), 1);
num_sui = size(find(subject_index == 3), 1);


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

%%
sv = svd(subj_coef_train);
scatter(1: 13, sv)

figure;
sv = svd(subj_coef_test);
scatter(1: 13, sv)

%%
clc
clc; clear;
load('files/cv_info.mat');
ctl_train_subs = cv_info.fold1.train.ctl;
sui_train_subs = cv_info.fold1.train.sui;

labels = [zeros(length(ctl_train_subs), 1); ones(length(sui_train_subs), 1)];
load('/hd2/research/EEG/code/output/ctl_sui_happy_fold1/subj_coef_happy_R2_ctl_sui_fold1.mat');
gscatter(subj_coef_train(:, 1), subj_coef_train(:, 2), labels);


