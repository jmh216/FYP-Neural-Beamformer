% TEST_PESQ2_MTLB demonstrates use of the PESQ2_MTLB function
%
%   See also PESQ2_MTLB
%
%   Author: Arkadiy Prodeus, email: aprodeus@gmail.com, July 2014

clear all; close all; clc;

% name of executable file for PESQ calculation
binary = 'pesq2.exe';
Audio_folder = 'C:\Users\74778\Desktop\PESQ\et05_bus_real_v2';

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

for j = 1:length(theFiles)/2
    clean = [];
    GEV = [];
    loop = 0;
    for i =  2*(j-1)+1 : 2*j      
        loop = loop+1;
        baseFileName = theFiles(i).name;
        fullFileName = fullfile(theFiles(i).folder, baseFileName);
        if loop == 2
            gev_name = baseFileName;
            gev =  v_readwav(fullFileName); 
        end
        if loop == 1
            clean_name = baseFileName;
            clean =  v_readwav(fullFileName);   
        end
    end
    
    expdir = pwd();
    pathaudio = Audio_folder;
    
    wb_GEV(j) = pesq2_mtlb( clean_name, gev_name, 16000, 'wb', binary, pathaudio);
    
    fprintf( 'Current Process  = %5.3f\n', j); 
end

fprintf( 'Current Process  = %5.3f\n', mean(wb_GEV,2)); 

