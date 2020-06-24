load('wb_dsb');
load('wb_mvdr');
load('wb_gev');

 wb_DSB_average = mean(wb_DSB,2);
 wb_MVDR_average = mean(wb_MVDR,2);
 wb_GEV_average = mean(wb_GEV,2);

fprintf('==================================\n'); 
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