clc; clear all; close all;
% https://uk.mathworks.com/help/stats/boxplot.html

load('DSB');
load('MVDR');
load('GEV');
wb_GEV_NOBAN = wb_GEV;
load('REF');
load('GEV+BAN')
wb_GEV_BAN = wb_GEV;
load('FW+BAN')
wb_FW_BAN = wb_GEV;
load('FW')
wb_FW = wb_GEV;


wb_DSB_average = mean(wb_DSB,2);
wb_MVDR_average = mean(wb_MVDR,2);
wb_REF_average = mean(wb_REF,2);
wb_GEV_average = mean(wb_GEV_NOBAN,2);
wb_GEVBAN_average = mean(wb_GEV_BAN,2);
wb_FWBAN_average = mean(wb_FW_BAN,2);
wb_FW_average = mean(wb_FW,2);


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
disp('WB-PESQ score for GEV:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_GEV_average );
fprintf('==================================\n'); 

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for GEV with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_GEVBAN_average );

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for FW with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_FWBAN_average );
fprintf('==================================\n'); 

    % display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for FW with BAN:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_FW_average );
fprintf('==================================\n');

PESQ_matrix_STR = [wb_REF', wb_DSB', wb_MVDR', wb_FW', wb_FW_BAN', wb_GEV_NOBAN',wb_GEV_BAN'];

figure
boxplot(PESQ_matrix_STR,'Labels',{'REF', 'DSB', 'MVDR', 'FF+GEV', 'FF+GEV+BAN', 'BLSTM+GEV', 'BLSTM+GEV+BAN'})
xlabel('Beamforming techniques')
ylabel('PESQ Score')

title('Box and whisker plot of PESQ for different beamformering techniques (For STR Real Data)')
