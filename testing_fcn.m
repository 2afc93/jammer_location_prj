clear all;
clc
close all
test_running = true;
test_uav_psi = 45*pi/180;
% test_jam_pos_array = 3000 + 6000*rand(2,10);
test_jam_pos_array = [3000 3000 3000 3000 3000 4500 4500 4500 4500 4500 6000 6000 6000 6000 6000 7500 7500 7500 7500 7500 9000 9000 9000 9000 9000;...
                      3000 4500 6000 7500 9000 3000 4500 6000 7500 9000 3000 4500 6000 7500 9000 3000 4500 6000 7500 9000 3000 4500 6000 7500 9000];
test_ranges = sqrt(sum(test_jam_pos_array.^2,1));
test_uav_pos = [900;900];
% test_filter = {'EKF', 'UKF', 'PF'};


P_t_min=50*(10^(-3));                                                %   [W] - Generally around 1mW
P_t_max=650*(10^(-3));                                              %   [W] - Generally 650mW
test_P_t_array=linspace(P_t_min, P_t_max, 11);

test_n_pos = size(test_jam_pos_array,2);
test_n_power = length(test_P_t_array);


test_error2 = nan(test_n_pos*test_n_power, 2401);
test_execTime = nan(test_n_power,test_n_pos);%zeros(length(test_P_t_array), length(test_jam_pos_array), 3);
test_Gtratio_error2 = nan(test_n_pos*test_n_power, 2401);
% test_filter = 'EKF';
% test_filter = 'UKF';

test_filter = 'EKF_aug';
test_Npart = 300;
for test_iter_Pt = 1:length(test_P_t_array)
    test_P_t = test_P_t_array(test_iter_Pt);
    for test_iter_pos = 1:length(test_jam_pos_array)
        test_jam_pos = test_jam_pos_array(:,test_iter_pos);
        run('Main_anisotropic');  
        test_range_evol = sqrt(sum(x_state.^2,1));
        test_error2(test_n_pos*(test_iter_Pt - 1) + test_iter_pos, :)           = (test_ranges(test_iter_pos)*ones(1,size(test_range_evol,2)) - test_range_evol).^2;
        test_execTime(test_iter_Pt, test_iter_pos)                              = nanmean(exec_time);
        test_Gtratio_error2(test_n_pos*(test_iter_Pt - 1) + test_iter_pos, :)   = (Gt_ratio' - system.Gt_ratio).^2;
    end
end

save('./results/anisotropic/EKF_aug_results', 'test_error2', 'test_execTime', 'test_Gtratio_error2');
disp('Analysis is completed')
clearvars -except test_*;

% test_error2 = nan(test_n_pos*test_n_power, 2401);
% test_execTime = nan(test_n_power,test_n_pos);
% test_Gtratio_error2 = nan(test_n_pos*test_n_power, 2401);
% 
% test_filter = 'EKF_aug';
% test_Npart = 500;
% for test_iter_Pt = 1:length(test_P_t_array)
%     test_P_t = test_P_t_array(test_iter_Pt);
%     for test_iter_pos = 1:length(test_jam_pos_array)
%         test_jam_pos = test_jam_pos_array(:,test_iter_pos);
%         
%         run('Main_anisotropic');  
%         test_range_evol = sqrt(sum(x_state.^2,1));
%         test_error2(test_n_pos*(test_iter_Pt - 1) + test_iter_pos, :) = (test_ranges(test_iter_pos)*ones(1,size(test_range_evol,2)) - test_range_evol).^2;
%         test_execTime(test_iter_Pt, test_iter_pos) = nanmean(exec_time);
%         test_Gtratio_error2(test_n_pos*(test_iter_Pt - 1) + test_iter_pos, :)   = (Gt_ratio - system.Gt_ratio)^2;
%     end
% end
% 
% save('./results/anisotropic/EKF_aug_results', 'test_error2', 'test_execTime', 'test_Gtratio_error2');
disp('Analysis is completed')
% clearvars -except test_*;