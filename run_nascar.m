clc; clear;
code_path = '/hd2/research/EEG/code';
addpath(genpath(code_path));

%%

load('/hd2/research/EEG/data/Darpa_Data_emo.mat');

%%

sad_data = Data.sad;
[n_channel, n_time, n_sub] = size(cur_data);
clearvars Data

num_control = size(find(subject_index == 1), 1);
num_depressed = size(find(subject_index == 2), 1);
num_suicidal = size(find(subject_index == 3), 1);

%% subject label
% 1: control, 2: depressed, 3: suicidal

n_per_class = 25;

sad_data2 = zeros(n_channel, n_time, n_per_class * 2);

cnt = 1;

for i = subject_index
    sad_data2(:, :, cnt) = sad_data(:, :, i);
    cnt = cnt + 1;
end
clearvars cnt

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
result = srscpd(sad_data2, R, option);

save('output/sad_R5_rst.mat');

