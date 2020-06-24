% TEST_PESQ2_MTLB demonstrates use of the PESQ2_MTLB function
%
%   See also PESQ2_MTLB
%
%   Author: Arkadiy Prodeus, email: aprodeus@gmail.com, July 2014

clear all; close all; clc;

% name of executable file for PESQ calculation
binary = 'pesq2.exe';
Audio_folder = 'C:\Users\74778\Desktop\PESQ\et05_bus_real';

if ~isfolder(Audio_folder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', Audio_folder);
    uiwait(warndlg(errorMessage));
    Audio_folder = uigetdir(); % Ask for a new one.
    if Audio_folder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(Audio_folder, '*.wav'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

for j = 1:length(theFiles)/5
    all = [];
    clean = [];
    DSB = [];
    MVDR = [];
    GEV = [];
    ref = [];
    loop = 0;
    for i =  5*(j-1)+1 : 5*j     
        loop = loop+1;
        baseFileName = theFiles(i).name;
        fullFileName = fullfile(theFiles(i).folder, baseFileName);
        if loop == 3
            gev_name = baseFileName;
            gev =  v_readwav(fullFileName); 
        end
        if loop == 4
            dsb_name = baseFileName;
            dsb =  v_readwav(fullFileName); 
        end
        if loop == 5
            mvdr_name = baseFileName;                
            mvdr =  v_readwav(fullFileName);      
        end
        if loop == 1
            clean_name = baseFileName;
            clean =  v_readwav(fullFileName);   
        end
        if loop == 2 
            ref_name = baseFileName;
            ref =  v_readwav(fullFileName);      
        end
    end
%     gev = gev(1:length(clean));
    
%     plot(gev);
%         plot(dsb);
%     hold on;
%     plot(dsb);
%     plot(mvdr);
%     plot(ref);  
%     plot(clean);  
%     plot(dsb-clean)
%     noise = (gev+10) - (clean+10);
%     plot(noise);
%     legend('1','2')
    
    expdir = pwd();
    pathaudio = Audio_folder;
    
    wb_REF(j) = pesq2_mtlb( clean_name, ref_name, 16000, 'wb', binary, pathaudio);
    wb_DSB(j) = pesq2_mtlb( clean_name, dsb_name, 16000, 'wb', binary, pathaudio);
    wb_MVDR(j) = pesq2_mtlb( clean_name, mvdr_name, 16000, 'wb', binary, pathaudio);
    wb_GEV(j) = pesq2_mtlb( clean_name, gev_name, 16000, 'wb', binary, pathaudio);
    
%     SNR_REF(j) = calSNR(clean,ref);
%     SNR_DSB(j) = calSNR(clean,dsb);
%     SNR_MVDR(j) = calSNR(clean,mvdr);
%     SNR_GEV(j) = calSNR(clean,gev);

    fprintf( 'Current Process  = %5.3f\n', j); 
end

wb_ref_average = mean(wb_REF,2);
wb_DSB_average = mean(wb_DSB,2);
wb_MVDR_average = mean(wb_MVDR,2);
wb_GEV_average = mean(wb_GEV,2);

fprintf('==================================\n'); 
% display results to screen
fprintf('====================\n'); 
disp('WB-PESQ score for reference channel:');
fprintf( 'WB MOS LQO  = %5.3f\n', wb_ref_average );

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


