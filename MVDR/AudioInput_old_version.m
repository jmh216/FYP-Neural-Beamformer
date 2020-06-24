%--------------------------------------------------------------------------
% Example script for DSB and MVDR response of audio input
% Author        : Jinming Hu
% Date          : 01-06-2020
% Version       : Version 2
%--------------------------------------------------------------------------
%% Parameter Settings
clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);

%% setup paths for mvdr
restoredefaultpath
addpath('submodules/sap-elobes-microphone-arrays');
addpath('submodules/sap-elobes-utilities');
addpath('submodules/sap-voicebox/voicebox');
addpath('lib');

%% beamformer parameters
% ensure empty struct
spec = struct();
% sample rate
spec.fs = 16000;            

% microphone array
nSensors = 6;
spacing = 0.15;
spec.ema_fcn = @()FreeField_20180202_01_ULA_configurable(nSensors,spacing); % function handle for elobes microphone array to use , create ULA position matrix 
                                                                                   
% Create cartesian unit vector for looking direction
steer_angle = 90; % 90 degree-[0 1 0] ; 0 degree-[1 0 0]
unit_source_position = [cosd(steer_angle) sind(steer_angle) 0]; 
spec.look_unit_vec = unit_source_position;  

%% MVDR weights
filtCoefs_mvdr = design_beamformer(spec);

%% Find RIR from looking direction
ema = spec.ema_fcn();
ema.setSoundSpeed(v_soundspeed);
ema.prepareData(spec.fs);
[target_az,target_inc,~] = mycart2sph(spec.look_unit_vec);
h = ema.getImpulseResponseForSrc(target_az,target_inc);

%% Take audio file as input

[speech,fs] = v_readwav('./submodules/sap-voicebox/voicebox_demo/data/as01b0.wav');
% [speech,fs] = v_readwav('female_speech.wav');

t = 0:1/fs:(length(speech)-1)/fs;
mic_speech = fftfilt(h,speech);

diffuse_noise = zeros(size(mic_speech));
noise_az_deg = (0:5:355).';
h_circ = ema.getImpulseResponseForSrc(deg2rad(noise_az_deg),deg2rad(90*ones(size(noise_az_deg))));
for idoa = 1:length(noise_az_deg)
    new_mono_noise = v_stdspectrum(11,'t',fs,length(speech));
    diffuse_noise = diffuse_noise + 1/8000 * fftfilt(squeeze(h_circ(:,:,idoa)),new_mono_noise);
end

mixed = mic_speech + diffuse_noise;

mixed_ref = mixed(:,ema.refChan);
target_out = sum(fftfilt(filtCoefs_mvdr,mic_speech),2);
bf_out = sum(fftfilt(filtCoefs_mvdr,mixed),2);

%% DSB
% setting
c = 343.2354656149216;                % Sound velocity (m/s)
n = 1600;               % Number of samples for RIR
beta = 0;               % Reverberation time (s)
Room = [8 8 4];			% Room dimensions [x y z] (m) 
Fs = 16000;             % Sample frequency (samples/s) 
T = 1/Fs;               % Sampling period (s)       
Mic_number = 6;         % Number of microphones
Mic_spacing = 0.15;     % Spacing of microphones (m)
steering_angle = 90;
Source_angle = 90;
Mic_Source_Distance = 1; % Distance between source and receiver
% ----- %
r = CreateReceiverPosition(Mic_number,Mic_spacing,Room(1),Room(2),Room(3)); % Generate receiver position [x y z] (m)
source_input_rir = mixed';

time_length = length(source_input_rir(1,:)); 
if mod(time_length,2) == 0 % N is even
    frequency_axis1 = (0:((time_length)/2))*Fs/(time_length); 
    frequency_axis2 = [frequency_axis1,fliplr(-frequency_axis1(2:end-1))];    
else % N is odd 
    frequency_axis1 = (0:((time_length-1)/2))*Fs/(time_length-1); 
    frequency_axis2 = [frequency_axis1,fliplr(-frequency_axis1(2:end))];
end

signal_input_freq = fft(source_input_rir.').';  % perform fft 
weights_DSB = FindWeightsDSB(r,Mic_spacing,frequency_axis2,steering_angle,c); % calculate weights of DSB
signal_output_DSB_freq = sum(weights_DSB.*signal_input_freq); % multiply the weights and add all the channels together
signal_output_DSB_time = real(ifft(signal_output_DSB_freq));

%% Plot

figure;
plot(t,mixed_ref);hold all;plot(t,bf_out);plot(t,target_out);
xlabel("Time/s"); ylabel("Amplitude");
legend('mixed reference','MVDR bf out','target out')
title('MVDR result');
hold off;

figure;
plot(t,mixed_ref);hold on;plot(t,signal_output_DSB_time);plot(t,target_out);
xlabel("Time/s"); ylabel("Amplitude");
legend('mixed reference','DSB out','target out');
title('DSB result');

%% Functions
function receiver_position = CreateReceiverPosition(Mic_number,Micro_Spacing,x,y,z)
    receiver_position = zeros(Mic_number,3); % A matrix with size nMic x 3, 3 for three axis x ,y, z
    Array_length = Micro_Spacing*(Mic_number-1);
    for index = 1:Mic_number
        receiver_position(index,:) = [(x-Array_length)/2+(index-1)*Micro_Spacing y/2 z/2]; % Microphones'position different in x axis
    end
end

function weights_DSB = FindWeightsDSB(Receiver_position, Mic_spacing, frequency_array, steering_angle, c)
    Mic_number = size(Receiver_position,1);
    frequency_upperbound = length(frequency_array);
    weights_DSB = zeros(Mic_number,frequency_upperbound);
    for Mic_index = 1:Mic_number  % for each microphone
        for frequency_index = 1:frequency_upperbound  % for each frequency
            weights_DSB(Mic_index,frequency_index) = 1/Mic_number*exp(-1*1i*frequency_array(frequency_index)*2*pi/c*Mic_spacing*(Mic_index-1)*cos(steering_angle*pi/180)); 
        end
    end
end

% function signal_output = ConvWithRIR(signal_input,RIR) 
%     Mic_number = size(RIR,1); 
%     sample_number = size(signal_input,2);
%     signal_output = zeros(Mic_number,sample_number);
%     for index = 1:Mic_number
%         signal_output(index,:) = filter(RIR(index,:),1,signal_input); 
%     end
% end
