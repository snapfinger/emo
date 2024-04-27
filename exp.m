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

U2 = result(1).U;
sub_coef = U2{3, 1};
labels = [zeros(15, 1); ones(15, 1)];
% gscatter(sub_coef(:, 1), sub_coef(:, 2), labels);
scatter(1: 30, sub_coef(:, 1));


