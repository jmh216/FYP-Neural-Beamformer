% TEST_PESQ2_MTLB demonstrates use of the PESQ2_MTLB function
%
%   See also PESQ2_MTLB
%
%   Author: Arkadiy Prodeus, email: aprodeus@gmail.com, July 2014

clear all; close all; clc;

% name of executable file for PESQ calculation
binary = 'pesq2.exe';
Audio_in_folder = 'C:\Users\74778\Desktop\Audio files clean\et05_bus_real';
Audio_out_folder_DSB = 'C:\Users\74778\Desktop\Audio files DSB\et05_bus_real';
Audio_out_folder_MVDR = 'C:\Users\74778\Desktop\Audio files MVDR\et05_bus_real';
Audio_out_folder_GEV = 'C:\Users\74778\Desktop\Audio files GEV\et05_bus_real';

if ~isfolder(Audio_in_folder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', Audio_in_folder);
    uiwait(warndlg(errorMessage));
    Audio_in_folder = uigetdir(); % Ask for a new one.
    if Audio_in_folder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(Audio_in_folder, '*.wav'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

filePattern_DSB = fullfile(Audio_out_folder_DSB, '*.wav'); % Change to whatever pattern you need.
theFiles_DSB = dir(filePattern_DSB);

filePattern_MVDR = fullfile(Audio_out_folder_MVDR, '*.wav'); % Change to whatever pattern you need.
theFiles_MVDR = dir(filePattern_MVDR);

filePattern_GEV = fullfile(Audio_out_folder_GEV, '*.wav'); % Change to whatever pattern you need.
theFiles_GEV = dir(filePattern_GEV);

for j = 1:length(theFiles)
    clean = [];
    DSB = [];
    MVDR = [];
    GEV = [];

    baseFileName = theFiles(j).name;
    fullFileName = fullfile(theFiles(j).folder, baseFileName);
    [clean, fs] =  v_readwav(fullFileName);    

    baseFileName_DSB = theFiles_DSB(j).name;
    fullFileName_DSB = fullfile(theFiles_DSB(j).folder, baseFileName_DSB);
    [DSB, fs] =  v_readwav(fullFileName_DSB);    

    baseFileName_MVDR = theFiles_MVDR(j).name;
    fullFileName_MVDR = fullfile(theFiles_MVDR(j).folder, baseFileName_MVDR);
    [MVDR, fs] =  v_readwav(fullFileName_MVDR);    

    baseFileName_GEV = theFiles_GEV(j).name;
    fullFileName_GEV = fullfile(theFiles_GEV(j).folder, baseFileName_GEV);
    [GEV, fs] =  v_readwav(fullFileName_GEV);    
    
    expdir = pwd();
    
    audiowrite([expdir, '\sounds\', baseFileName], clean, fs)
    audiowrite([expdir, '\sounds\', baseFileName_DSB], DSB, fs)
    audiowrite([expdir, '\sounds\', baseFileName_MVDR], MVDR, fs)
    audiowrite([expdir, '\sounds\', baseFileName_GEV], GEV, fs)
    
    %% PESQ
         
    pathaudio = 'sounds';
    wb_DSB(j) = pesq2_mtlb( baseFileName, baseFileName_DSB, 16000, 'wb', binary,pathaudio);
    wb_MVDR(j) = pesq2_mtlb( baseFileName, baseFileName_MVDR, 16000, 'wb', binary,pathaudio);
    wb_GEV(j) = pesq2_mtlb( baseFileName, baseFileName_GEV, 16000, 'wb', binary,pathaudio);

end

wb_DSB_average = mean2(wb_DSB);
wb_MVDR_average = mean2(wb_MVDR);
wb_GEV_average = mean2(wb_GEV);

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


