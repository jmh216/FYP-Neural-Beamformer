%--------------------------------------------------------------------------
% Example script for DSB response of audio input
% Author        : Jinming Hu
% Date          : 05-06-2020
% Version       : Final
%--------------------------------------------------------------------------
%% Parameter Settings
clc;
% close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);

c = 340;                % Sound velocity (m/s)
n = 1024;               % Number of samples for RIR
beta = 0;               % Reverberation time (s)
Room = [8 8 4];			% Room dimensions [x y z] (m) 

Fs = 16000;             % Sample frequency (samples/s) 
T = 1/Fs;               % Sampling period (s)       

Mic_number = 6;         % Number of microphones
Mic_spacing = 0.15;     % Spacing of microphones (m)
Mic_Source_Distance = 2; % Distance between source and receiver

steering_angle = 0;
source_angle = 0;
noise_angle = 0;

% read audio input 
filename_input = 'female_speech.wav';
[audio_in, fs] = audioread(filename_input); %sampled data, y, and a sample rate for that data, Fs=16000

filename_inter = 'clean.wav';
[audio_inter] = audioread(filename_inter,[137549, 275097]); %sampled data, y, and a sample rate for that data, Fs=16000

load('audio_in_rir.mat');

load('0-0-0.mat');
signal_output_DSB_time_0 = signal_output_DSB_time;

load('0-0-90.mat');
signal_output_DSB_time_90 = signal_output_DSB_time;

t = 0:1/fs:(length(audio_in)-1)/fs;

figure;
subplot(221);plot(t,audio_in);
xlabel("Time/s"); ylabel("Amplitude");
ylim([-1,1]);
title("Original target signal waveform");

subplot(222);plot(t,audio_inter);
xlabel("Time/s"); ylabel("Amplitude");
ylim([-1,1]);
title("Original interference signal waveform");

subplot(223);plot(t,signal_output_DSB_time_0);
hold on;
plot(t,audio_in_rir(1,:));
legend('DSB Audio out','Ref Audio in');
xlabel("Time/s"); ylabel("Amplitude");
ylim([-0.05,0.05]);
title("Steer at 0 - Target at 0 - Interference at 0 degree");
hold off;

subplot(224);plot(t,signal_output_DSB_time_90);
hold on;
plot(t,audio_in_rir(1,:));
legend('DSB Audio out','Ref Audio in');
xlabel("Time/s"); ylabel("Amplitude");
ylim([-0.05,0.05]);
title("Steer at 0 - Target at 0 - Interference at 90 degree");
hold off;

figure;
s1 = spectrogram(signal_output_DSB_time_0);
spectrogram(signal_output_DSB_time_0,'yaxis')

% figure;
% s2 = spectrogram(signal_output_DSB_time_90);
% spectrogram(signal_output_DSB_time_90,'yaxis')

figure;
v_spgrambw(signal_output_DSB_time_0,fs,'pJcw') ;
