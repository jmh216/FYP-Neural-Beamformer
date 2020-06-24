clc; clear all; close all;
% https://uk.mathworks.com/help/stats/boxplot.html
set(0,'DefaultAxesFontSize',13)
set(0,'DefaultLineLineWidth',1);

load('BUS');
load('CAF');
load('PED');
load('STR');

PEDQ_ALL = [PESQ_matrix__bus([1:100],:); PESQ_matrix_caf;  PESQ_matrix_ped;  PESQ_matrix_STR];

wb_REF_average = mean(PEDQ_ALL(:,1),1);
wb_DSB_average = mean(PEDQ_ALL(:,2),1);
wb_MVDR_average = mean(PEDQ_ALL(:,3),1);
wb_FW_average = mean(PEDQ_ALL(:,4),1);
wb_FWBAN_average = mean(PEDQ_ALL(:,5),1);
wb_GEV_average = mean(PEDQ_ALL(:,6),1);
wb_GEVBAN_average = mean(PEDQ_ALL(:,7),1);

fprintf('==================================\n'); 
% display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for REF:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_REF_average );

% display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for DSB:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_DSB_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for MVDR:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_MVDR_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for FW with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_FW_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for FW with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_FWBAN_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for GEV:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_GEV_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for GEV with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_GEVBAN_average );
fprintf('==================================\n'); 

figure
boxplot(PEDQ_ALL,'Labels',{'REF', 'DSB', 'MVDR', 'FF+GEV', 'FF+GEV+BAN', 'BLSTM+GEV', 'BLSTM+GEV+BAN'})
xlabel('Beamforming techniques')
ylabel('PESQ Score')

title('Box and whisker plot of PESQ for different beamformering techniques (Average)')
