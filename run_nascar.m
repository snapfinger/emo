clc; clear;
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));
load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

%% 
num_ctl = size(find(subject_index == 1), 1);
num_dep = size(find(subject_index == 2), 1);
num_sui = size(find(subject_index == 3), 1);
%     
load('files/subject_type_list.mat');
load('files/cv_info.mat');

%%

% construct the training data
% subject label:
% 1: control, 2: depressed, 3: suicidal

% NASCAR rank
R = 20;
    
for fold = 1: 5
% for fold = 1
    
    fold_name = ['fold' num2str(fold)]
    
    for task = ["happy", "sad", "neutral", "suicide"]
%     for task = ["happy"]

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
        
        cur_fold_info = cv_info.(fold_name);
        
        train_subs_ctl = cur_fold_info.train.ctl;
        train_subs_sui = cur_fold_info.train.sui;
        n_train = length(train_subs_ctl) + length(train_subs_sui);
        
        cur_train_data2 = zeros(n_channel, n_time, n_train);
        ctl_train_grand_sum = zeros(n_channel, n_time);
        
        cnt_train = 0;
        for sub = train_subs_ctl
            cnt_train = cnt_train + 1;
            cur_train_data2(:, :, cnt_train) = cur_data(:, :, sub);
            ctl_train_grand_sum = ctl_train_grand_sum + cur_data(:, :, sub);
        end

        ctl_train_grand_avg = ctl_train_grand_sum / length(train_subs_ctl);
        clearvars ctl_train_grand_sum

        for sub = train_subs_sui
            cnt_train = cnt_train + 1;
            cur_train_data2(:, :, cnt_train) = cur_data(:, :, sub); 
        end

        % subtract control train grand average (#channel * #time) from each subject's data
        for i = 1: n_train
            cur_train_data2(:, :, i) = cur_train_data2(:, :, i) - ctl_train_grand_avg;
        end

        option = srscpd('opt');
        option.nonnegative = [0 0 0];
        option.rankOneOptALS.useTFOCS = false; % if use TFOCS or not to solve rank 1 problem
        option.optAlg.normRegParam = 0.001; % regularization parameter on the norm of each mode
        option.optAlg.optSolver.printItv = 50; % print interval
        option.optAlg.optSolver.learningRate = 1e-3; % Nadam learning rate
        option.optAlg.optSolver.maxNumIter = 2000; % max number of iterations in Nadam
        option.optAlg.cacheMTS = true; % turn this on if you have enough memory to speed up
        option.logFile = 'log_file';
        result = srscpd(cur_train_data2, R, option);

        rst_folder = ['output/ctl_sui_' str2char(task) '_fold' num2str(fold) '/'];
        if ~exist(rst_folder, 'dir')
            mkdir(rst_folder)
        end
        
        save([rst_folder 'nascar_' str2char(task) '_R' num2str(R) ...
                '_ctl_sui_fold' num2str(fold) '.mat'], 'result');

        % fitting for test subjects
        for cur_R = 2: R 

            basis = zeros(n_channel * n_time, cur_R);
            U_rst = result(cur_R).U;
            nascar_spatial = U_rst{1};
            nascar_temporal = U_rst{2};
            subj_coef_train = U_rst{3};

            for i = 1: cur_R
                cur_temporal_spatio = nascar_spatial(:, i) .* nascar_temporal(:, i)';
                basis(:, i) = cur_temporal_spatio(:);
            end

            test_subs_ctl = cur_fold_info.test.ctl;
            test_subs_sui = cur_fold_info.test.sui;
            n_test = length(test_subs_ctl) + length(test_subs_sui);
        
            cur_test_data2 = zeros(n_channel, n_time, n_test);

            cnt_test = 0;
            for sub = test_subs_ctl
                cnt_test = cnt_test + 1;
                cur_test_data2(:, :, cnt_test) = ...
                    cur_data(:, :, sub) - ctl_train_grand_avg;
            end

            for sub = test_subs_sui
                cnt_test = cnt_test + 1;
                cur_test_data2(:, :, cnt_test) = ...
                    cur_data(:, :, sub) - ctl_train_grand_avg;
            end

            subj_coef_test = zeros(cur_R, n_test);

            for i = 1: n_test
                cur_eeg = cur_test_data2(:, :, i);
                cur_eeg = cur_eeg(:);
                subj_coef_test(:, i) = basis \ cur_eeg;
            end

            subj_coef_test = subj_coef_test';
            subj_coef_test = subj_coef_test ./ vecnorm(subj_coef_test);

            save([rst_folder 'subj_coef_' str2char(task) '_R' num2str(cur_R) ...
                    '_ctl_sui_fold' num2str(fold) '.mat' ], ...
                    'subj_coef_train', 'subj_coef_test')
        end
    end 
end






