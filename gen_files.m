%%
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

%%
clc; 
load('files/subject_type_list.mat');

% fold 1: 
% test: 
%    ctl: [1 - 7], sui: [1 - 5], dep: [1 - 4]

% fold 2
% test:
%   ctl: [8 - 14], sui: [6 - 10], dep: [5 - 8]

% fold 3
% test
%   ctl: [15 - 21], sui: [11 - 15], dep: [9 - 13]

% fold 4
% test
%   ctl: [22 - 28], sui: [16 - 20], dep: [14 - 18]

% fold 5
% test
%   ctl: [29 - 35], sui: [21 - 26], dep: [19 - 23]

num_ctl = length(ctl_id_list);
num_dep = length(dep_id_list);
num_sui = length(sui_id_list);


cv_info = struct();

for fold = 1: 5
    
    fold

    if fold == 1
        test_ctl_id = 1: 7;
        test_sui_id = 1: 5;
        test_dep_id = 1: 4;

    elseif fold == 2
        test_ctl_id = 8: 14;
        test_sui_id = 6: 10;
        test_dep_id = 5: 8;
    elseif fold == 3
        test_ctl_id = 15: 21;
        test_sui_id = 11: 15;
        test_dep_id = 9: 13;
    elseif fold == 4
        test_ctl_id = 22: 28;
        test_sui_id = 16: 20;
        test_dep_id = 14: 18;
    elseif fold == 5
        test_ctl_id = 29: 35;
        test_sui_id = 21: 26;
        test_dep_id = 19: 23;
    end

    train_ctl_id = setdiff(1: num_ctl, test_ctl_id);
    train_sui_id = setdiff(1: num_sui, test_sui_id);
    train_dep_id = setdiff(1: num_dep, test_dep_id);

    test_ctl_subs = ctl_id_list(test_ctl_id);
    train_ctl_subs = ctl_id_list(train_ctl_id);

    test_sui_subs = sui_id_list(test_sui_id);
    train_sui_subs = sui_id_list(train_sui_id);

    test_dep_subs = dep_id_list(test_dep_id);
    train_dep_subs = dep_id_list(train_dep_id);

    train_list = struct('ctl', train_ctl_subs, 'sui', train_sui_subs, 'dep', train_dep_subs);
    test_list = struct('ctl', test_ctl_subs, 'sui', test_sui_subs, 'dep', test_dep_subs);
    
    cur_fold = struct('train', train_list, 'test', test_list);
    cur_fold_name = ['fold' num2str(fold)];
    cv_info.(cur_fold_name) = cur_fold;
    
end


save('files/cv_info.mat', 'cv_info');






