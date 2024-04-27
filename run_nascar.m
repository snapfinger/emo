clc; clear;
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));
load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

%%

task = 'happy'

if strcmp(task, 'happy')
    cur_data = Data.happy;
elseif strcmp(task, 'sad')
    cur_data = Data.sad;
end

%%
[n_channel, n_time, n_sub] = size(cur_data);
clearvars Data

num_ctl = size(find(subject_index == 1), 1);
num_dep = size(find(subject_index == 2), 1);
num_sui = size(find(subject_index == 3), 1);

%% form the training data
% choose 25 control, 25 suicidal
% subject label
% 1: control, 2: depressed, 3: suicidal


clc

n_per_class = 15;
cur_data2 = zeros(n_channel, n_time, n_per_class * 2);
ctl_grand_avg = zeros(n_channel, n_time);

cnt_ctl = 1;
for i = 1: length(subject_index)
    fprintf("i=%d, cnt_ctl=%d\n", i, cnt_ctl);
    
    if cnt_ctl == n_per_class + 1
        break;
    end
    
    if subject_index(i) == 1
        cur_data2(:, :, cnt_ctl) = cur_data(:, :, i);
        ctl_grand_avg = ctl_grand_avg + cur_data(:, :, i);
        cnt_ctl = cnt_ctl + 1;
    end    
end

ctl_train_grand_avg = ctl_grand_avg / n_per_class;

%%
cnt_dep = 1;
for i = 1: length(subject_index)
    fprintf("i=%d, cnt_dep=%d\n", i, cnt_dep);
    
    if cnt_dep == n_per_class + 1
        break;
    end
    
    if subject_index(i) == 3
        cur_data2(:, :, n_per_class + cnt_dep) = cur_data(:, :, i);
        cnt_dep = cnt_dep + 1;
    end 
end

clearvars cnt

%% subtract control train grand average (#channel * #time) from each subject's data

for i = 1: n_per_class * 2
    cur_data2(:, :, i) = cur_data2(:, :, i) - ctl_train_grand_avg;
end

%%

R = 3;

option = srscpd('opt');
option.nonnegative = [0 0 1];
option.rankOneOptALS.useTFOCS = false; % if use TFOCS or not to solve rank 1 problem
option.optAlg.normRegParam = 0.001; % regularization parameter on the norm of each mode
option.optAlg.optSolver.printItv = 10; % print interval
option.optAlg.optSolver.learningRate = 1e-3; % Nadam learning rate
option.optAlg.optSolver.maxNumIter = 2000; % max number of iterations in Nadam
option.optAlg.cacheMTS = true; % turn this on if you have enough memory to speed up
% option.saveToFile = 'file_to_save_result';
% save logs (what printed in the command window) to a file
option.logFile = 'log_file';
result = srscpd(cur_data2, R, option);

%%
save(['output/' task '_R' num2str(R) '_rst.mat'], 'result');




