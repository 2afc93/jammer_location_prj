clear all; close all;
RMSE_lim = 40;
error_reg = [0.1 1.0 3.0 5.0]';
k2end = 200;
%% UKF results
load('UKF_results.mat')
% UKF_sse = sqrt(test_error2(:,end));
% UKF_ferror_idx = UKF_sse<RMSE_lim;
% UKF_ferror = test_error2(UKF_ferror_idx, :);
% UKF_rmse = sqrt(mean(UKF_ferror,1));

UKF_sse = sqrt(mean(test_error2(:,end-k2end:end),2));
UKF_ferror_idx = UKF_sse<RMSE_lim;
UKF_ferror = test_error2(UKF_ferror_idx, :);
UKF_rmse = sqrt(mean(UKF_ferror,1));

UKF_perc(1) = sum(UKF_sse < error_reg(1));
UKF_perc(2) = sum((UKF_sse > error_reg(1)).*(UKF_sse < error_reg(2)));
UKF_perc(3) = sum((UKF_sse > error_reg(2)).*(UKF_sse < error_reg(3)));
UKF_perc(4) = sum((UKF_sse > error_reg(3)).*(UKF_sse < error_reg(4)));
% UKF_perc(5) = sum((UKF_sse > error_reg(4)).*(UKF_sse < error_reg(5)));
% UKF_perc(6) = sum((UKF_sse > error_reg(5)).*(UKF_sse < error_reg(6)));
UKF_perc(5) = sum((UKF_sse > error_reg(4)));
UKF_perc = UKF_perc/sum(UKF_perc)*100;

UKF_meanTime = mean(mean(test_execTime));

% figure(1)
% semilogy((1:1801)*0.5, UKF_rmse, 'Linewidth', 2.0)
% hold on;
% % semilogy((1:1801)*0.5, UKF_rmse+UKF_td, '--','Linewidth', 2.0)
% % semilogy((1:1801)*0.5, UKF_rmse-UKF_td, '--','Linewidth', 2.0)
% grid on;
% xlabel('Time [sec]', 'Fontsize', 14);
% ylabel('RMSE [m]', 'Fontsize', 14);
% xlim([0 1801*0.5]);
clear test_error2 test_execTime

%% EKF results
load('EKF_results.mat')
% EKF_sse = sqrt(test_error2(:,end));
% EKF_ferror_idx = EKF_sse<RMSE_lim;
% EKF_ferror = test_error2(EKF_ferror_idx, :);
% EKF_rmse = sqrt(mean(EKF_ferror,1));

EKF_sse = sqrt(mean(test_error2(:,end-k2end:end),2));
EKF_ferror_idx = EKF_sse<RMSE_lim;
EKF_ferror = test_error2(EKF_ferror_idx, :);
EKF_rmse = sqrt(mean(EKF_ferror,1));

EKF_perc(1) = sum(EKF_sse < error_reg(1));
EKF_perc(2) = sum((EKF_sse > error_reg(1)).*(EKF_sse < error_reg(2)));
EKF_perc(3) = sum((EKF_sse > error_reg(2)).*(EKF_sse < error_reg(3)));
EKF_perc(4) = sum((EKF_sse > error_reg(3)).*(EKF_sse < error_reg(4)));
% EKF_perc(5) = sum((EKF_sse > error_reg(4)).*(EKF_sse < error_reg(5)));
% EKF_perc(6) = sum((EKF_sse > error_reg(5)).*(EKF_sse < error_reg(6)));
EKF_perc(5) = sum((EKF_sse > error_reg(4)));
EKF_perc = EKF_perc/sum(EKF_perc)*100;

EKF_meanTime = mean(mean(test_execTime));

% figure(2)
% semilogy((1:1801)*0.5, EKF_rmse, 'Linewidth', 2.0)
% hold on; 
% grid on;
% xlabel('Time [sec]', 'Fontsize', 14);
% ylabel('RMSE [m]', 'Fontsize', 14);
% xlim([0 1801*0.5]);
clear test_error2 test_execTime

%% Particle Filter reults N=150
load('PF150_results.mat')
% PF150_sse = sqrt(test_error2(:,end));
% PF150_ferror_idx = PF150_sse<RMSE_lim;
% PF150_ferror = test_error2(PF150_ferror_idx, :);
% PF150_rmse = sqrt(nanmean(PF150_ferror,1));

PF150_sse = sqrt(mean(test_error2(:,end-k2end:end),2));
PF150_ferror_idx = PF150_sse<RMSE_lim;
PF150_ferror = test_error2(PF150_ferror_idx, :);
PF150_rmse = sqrt(nanmean(PF150_ferror,1));

PF150_perc(1) = sum(PF150_sse < error_reg(1));
PF150_perc(2) = sum((PF150_sse > error_reg(1)).*(PF150_sse < error_reg(2)));
PF150_perc(3) = sum((PF150_sse > error_reg(2)).*(PF150_sse < error_reg(3)));
PF150_perc(4) = sum((PF150_sse > error_reg(3)).*(PF150_sse < error_reg(4)));
% PF150_perc(5) = sum((PF150_sse > error_reg(4)).*(PF150_sse < error_reg(5)));
% PF150_perc(6) = sum((PF150_sse > error_reg(5)).*(PF150_sse < error_reg(6)));
PF150_perc(5) = sum((PF150_sse > error_reg(4)));

PF150_perc = PF150_perc/sum(PF150_perc)*100;

PF150_meanTime = mean(mean(test_execTime));

% figure(3)
% semilogy((1:1801)*0.5, PF150_rmse, 'Linewidth', 2.0)
% hold on; 
% grid on;
% xlabel('Time [sec]', 'Fontsize', 14);
% ylabel('RMSE [m]', 'Fontsize', 14);
% xlim([0 1801*0.5]);
clear test_error2 test_execTime

%% Particle Filter results N=300
load('PF300_results.mat')
% PF300_sse = sqrt(test_error2(:,end));
% PF300_ferror_idx = PF300_sse<RMSE_lim;
% PF300_ferror = test_error2(PF300_ferror_idx, :);
% PF300_rmse = sqrt(nanmean(PF300_ferror,1));

PF300_sse = sqrt(mean(test_error2(:,end-k2end:end),2));
PF300_ferror_idx = PF300_sse<RMSE_lim;
PF300_ferror = test_error2(PF300_ferror_idx, :);
PF300_rmse = sqrt(nanmean(PF300_ferror,1));

PF300_perc(1) = sum(PF300_sse < error_reg(1));
PF300_perc(2) = sum((PF300_sse > error_reg(1)).*(PF300_sse < error_reg(2)));
PF300_perc(3) = sum((PF300_sse > error_reg(2)).*(PF300_sse < error_reg(3)));
PF300_perc(4) = sum((PF300_sse > error_reg(3)).*(PF300_sse < error_reg(4)));
% PF300_perc(5) = sum((PF300_sse > error_reg(4)).*(PF300_sse < error_reg(5)));
% PF300_perc(6) = sum((PF300_sse > error_reg(5)).*(PF300_sse < error_reg(6)));
PF300_perc(5) = sum((PF300_sse > error_reg(4)));
PF300_perc = PF300_perc/sum(PF300_perc)*100;

PF300_meanTime = mean(mean(test_execTime));

% figure(4)
% semilogy((1:1801)*0.5, PF300_rmse, 'Linewidth', 2.0)
% hold on; 
% grid on;
% xlabel('Time [sec]', 'Fontsize', 14);
% ylabel('RMSE [m]', 'Fontsize', 14);
% xlim([0 1801*0.5]);
clear test_error2 test_execTime

%% Particle Filter results N=500
load('PF500_results.mat')
% PF500_sse = sqrt(test_error2(:,end));
% PF500_ferror_idx = PF500_sse<RMSE_lim;
% PF500_ferror = test_error2(PF500_ferror_idx, :);
% PF500_rmse = sqrt(nanmean(PF500_ferror,1));

PF500_sse = sqrt(mean(test_error2(:,end-k2end:end),2));
PF500_ferror_idx = PF500_sse<RMSE_lim;
PF500_ferror = test_error2(PF500_ferror_idx, :);
PF500_rmse = sqrt(nanmean(PF500_ferror,1));

PF500_perc(1) = sum(PF500_sse < error_reg(1));
PF500_perc(2) = sum((PF500_sse > error_reg(1)).*(PF500_sse < error_reg(2)));
PF500_perc(3) = sum((PF500_sse > error_reg(2)).*(PF500_sse < error_reg(3)));
PF500_perc(4) = sum((PF500_sse > error_reg(3)).*(PF500_sse < error_reg(4)));
% PF500_perc(5) = sum((PF500_sse > error_reg(4)).*(PF500_sse < error_reg(5)));
% PF500_perc(6) = sum((PF500_sse > error_reg(5)).*(PF500_sse < error_reg(6)));
PF500_perc(5) = sum((PF500_sse > error_reg(4)));
PF500_perc = PF500_perc/sum(PF500_perc)*100;

PF500_meanTime = mean(mean(test_execTime));

% figure(5)
% semilogy((1:1801)*0.5, PF500_rmse, 'Linewidth', 2.0)
% hold on; 
% grid on;
% xlabel('Time [sec]', 'Fontsize', 14);
% ylabel('RMSE [m]', 'Fontsize', 14);
% xlim([0 1801*0.5]);
clear test_error2 test_execTime

%% RMSE Comparison
figure(1)
legend_rmse{1} = 'Extended Kalman Filter';
legend_rmse{2} = 'Unscented Kalman Filter';
legend_rmse{3} = 'Particle Filter (N = 150)';
legend_rmse{4} = 'Particle Filter (N = 300)';
legend_rmse{5} = 'Particle Filter (N = 500)';
semilogy((1:1801)*0.5/60, EKF_rmse, 'Linewidth', 2.0)
hold on;
semilogy((1:1801)*0.5/60, UKF_rmse, 'Linewidth', 2.0)
hold on;
semilogy((1:1801)*0.5/60, PF150_rmse, 'Linewidth', 2.0)
hold on;
semilogy((1:1801)*0.5/60, PF300_rmse, 'Linewidth', 2.0)
hold on;
semilogy((1:1801)*0.5/60, PF500_rmse, 'Linewidth', 2.0)
grid on;
xlabel('Time [min]', 'Fontsize', 14);
ylabel('RMSE [m]', 'Fontsize', 14);
set(gca, 'FontSize', 14.0);
xlim([0 1801*0.5/60]);
legend(legend_rmse, 'Location', 'NorthEast', 'Fontsize', 14);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 40 20];

iptsetpref('ImshowBorder','tight');
print('-r300', '-depsc', '../../report/figures/rmse_isot')

%% Execution Time Comparison
figure(2)
bar(1:5,[EKF_meanTime, UKF_meanTime, PF150_meanTime, PF300_meanTime, PF500_meanTime]*1e3, 'FaceColor',[0 102 204]/255,'EdgeColor',[0 0 255]/255,'LineWidth',1.5)
hold on;
set(gca,'YScale','log','XTickLabel',legend_rmse)
set(gca, 'FontSize', 14.0);
grid on;
ylabel('Execution Time [msec]')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 40 20];

iptsetpref('ImshowBorder','tight');
print('-r300', '-depsc', '../../report/figures/exec_time_isot')

%% Success Comparison
rmse_legend{1} = ['RMSE < ', num2str(error_reg(1)), ' m'];
rmse_legend{2} = [num2str(error_reg(1)), ' m', ' < RMSE < ', num2str(error_reg(2)), ' m'];
rmse_legend{3} = [num2str(error_reg(2)), ' m', ' < RMSE < ', num2str(error_reg(3)), ' m'];
rmse_legend{4} = [num2str(error_reg(3)), ' m', ' < RMSE < ', num2str(error_reg(4)), ' m'];
% rmse_legend{5} = [num2str(error_reg(4)), ' m', ' < RMSE < ', num2str(error_reg(5)), ' m'];
% rmse_legend{6} = [num2str(error_reg(5)), ' m', ' < RMSE < ', num2str(error_reg(6)), ' m'];
rmse_legend{5} = [num2str(error_reg(4)), ' m', ' < RMSE'];


figure(3)
suc_matrix = [EKF_perc; UKF_perc; PF150_perc; PF300_perc; PF500_perc];
bar(suc_matrix, 'stacked')
set(gca,'FontSize',14.0, 'XTickLabel',legend_rmse)
hold on;
% set(gca,'YScale','log','XTickLabel',legend_rmse, 'XLabel')
grid on;
ylabel('%')
ylim([0 100])
legend(rmse_legend, 'Location','northoutside', 'FontSize', 14)

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 40 20];

iptsetpref('ImshowBorder','tight');
print('-r300', '-depsc', '../../report/figures/error_dist_isot')
% close all;