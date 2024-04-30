clc; clear;
addpath(genpath('/hd2/sw1/eeglab2024.0'));
load('/hd2/research/EEG/data/chanlocs.mat');

%%
clc
task = 'happy'; 
fold = 1;
R = 20;

load(['/hd2/research/EEG/code/output/ctl_sui_happy_fold1/nascar_' task ...
        '_R' num2str(R) '_ctl_sui_fold' num2str(fold) '.mat']);
    
U = result(R).U;
spatial = U{1, 1};
temporal = U{2, 1};
subject = U{3, 1};

figs_folder = '/hd2/research/EEG/code/figs/';

for i = 1: R
    
    f = figure;
    subplot(2,2,1);
    topoplot(spatial(:, i), chanlocs, 'maplimits', 'maxmin', 'headrad', 'rim');
    colorbar;
    title('Spatial mode')

    subplot(2,2,2);
    bar(subject(:, i));
    xlabel('subject');
    ylabel('participation level');
    title('Subject mode')

    subplot(2,2,[3,4]);
    plot(1: size(temporal, 1), temporal(:, i), 'LineWidth', 2);
    xlabel('#Time');
    title('Temporal mode')
    
    saveas(f, ['figs/' task '_fold' num2str(fold) '_net' num2str(i) '.png']);
    
    close(f);
end