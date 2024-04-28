clc; clear;

%%
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));
load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

num_ctl = size(find(subject_index == 1), 1);
num_dep = size(find(subject_index == 2), 1);
num_sui = size(find(subject_index == 3), 1);
    
%%

load('/hd2/research/EEG/data/subject_type_list.mat');

%%

% construct the training data: 15 control, 15 suicidal
% subject label:
% 1: control, 2: depressed, 3: suicidal

n_ctl_train = 28;
n_sui_train = 20;
n_dep_train = 0;
n_train = n_ctl_train + n_sui_train + n_dep_train;
    
    
% for task = ["happy", "sad", "neutral", "suicide"]
for task = ["happy"]
    
    task

    if strcmp(task, 'happy')
        cur_data = Data.happy;
    elseif strcmp(task, 'sad')
        cur_data = Data.sad;
    elseif strcmp(task, 'neutral')
        cur_data = Data.neutral;
    elseif strcmp(task, 'suicide')
        cur_data = Data.suicide;
    end

    
    [n_channel, n_time, n_sub] = size(cur_data);

   
    cur_data2 = zeros(n_channel, n_time, n_train);
    ctl_train_grand_sum = zeros(n_channel, n_time);

    for i = 1: n_ctl_train
        cur_data2(:, :, i) = cur_data(:, :, ctl_id_list(i));
        ctl_train_grand_sum = ctl_train_grand_sum + cur_data(:, :, ctl_id_list(i));
    end

    ctl_train_grand_avg = ctl_train_grand_sum / n_ctl_train;
    clearvars ctl_train_grand_sum

    
    % suicidal (3)
    for i = 1: n_sui_train
        cur_data2(:, :, n_ctl_train + i) = cur_data(:, :, sui_id_list(i));
       
    end

    clearvars cnt

    % subtract control train grand average (#channel * #time) from each subject's data
    for i = 1: n_train
        cur_data2(:, :, i) = cur_data2(:, :, i) - ctl_train_grand_avg;
    end

    % run NASCAR
    R = 3;

    option = srscpd('opt');
    option.nonnegative = [0 0 1];
    option.rankOneOptALS.useTFOCS = false; % if use TFOCS or not to solve rank 1 problem
    option.optAlg.normRegParam = 0.001; % regularization parameter on the norm of each mode
    option.optAlg.optSolver.printItv = 10; % print interval
    option.optAlg.optSolver.learningRate = 1e-3; % Nadam learning rate
    option.optAlg.optSolver.maxNumIter = 2000; % max number of iterations in Nadam
    option.optAlg.cacheMTS = true; % turn this on if you have enough memory to speed up
    option.logFile = 'log_file';
    result = srscpd(cur_data2, R, option);

    
    save(['output/' str2char(task) '_R' num2str(R) ...
            '_ctl' num2str(n_ctl_train) ...
            '_sui' num2str(n_sui_train) ...
            '_dep' num2str(n_dep_train) '_rst_subNone.mat'], ...
            'result');
end

%% fitting for test subjects
clc

n_ctl_test = 7;
n_sui_test = 6;
n_dep_test = 0;
n_test = n_ctl_test + n_sui_test + n_dep_test;
cur_test_data2 = zeros(n_channel, n_time, n_test);
cnt_test_ctl = 1;

for i = n_ctl_train: length(ctl_id_list)
    cur_test_data2(:, :, cnt_test_ctl) = cur_data(:, :, i);
end

basis = zeros(n_channel * n_time, R);

U2 = result(R).U;
nascar_spatial = U2{1};
nascar_temporal = U2{2};

for i = 1: R
    cur_temporal_spatio = nascar_spatial(:, i) .* nascar_temporal(:, i)';
    basis(:, i) = cur_temporal_spatio(:);
end

subj_coef = zeros(R, 40);

for i = 1: 40
    cur_eeg = cur_data2(:, :, i);
    cur_eeg = cur_eeg(:);
    subj_coef(:, i) = basis \ cur_eeg;
end

subj_coef = subj_coef';
subj_coef = subj_coef ./ vecnorm(subj_coef);

%%
% clc
% subj_coef_gt = U2{3, 1};
% norm(subj_coef - subj_coef_gt) / norm(subj_coef_gt)



