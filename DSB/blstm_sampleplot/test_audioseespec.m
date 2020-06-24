clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLineWidth',1);

filename_output1 = '1.wav';
    [audio_in1, fs] = audioread(filename_output1); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output2 = '2.wav';
    [audio_in2, fs] = audioread(filename_output2); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output5 = 'clean.wav';
    [audio_in5, fs] = audioread(filename_output5); %sampled data, y, and a sample rate for that data, Fs=16000
    
t = 0:1/fs:(length(audio_in1)-1)/fs;

figure;
subplot(221);
plot(t,audio_in1(1:length(t),1));
hold on;
plot(t,audio_in5(1:length(t),1))
xlabel('Time/s');ylabel('Amplitude');
legend('REF','Clean');
hold off;

subplot(222);
v_spgrambw(audio_in1,fs,'pJcw') ;
title('REF')

subplot(223);
plot(t,audio_in2(1:length(t),1));
hold on;
plot(t,audio_in5(1:length(t),1))
xlabel('Time/s');ylabel('Amplitude');
legend('BLSTM+GEV+BAN','Clean');
hold off;

subplot(224);
v_spgrambw(audio_in2,fs,'pJcw') ;
title('BLSTM+GEV+BAN')

% figure;
% title('Spectrogram Comparison');

figure;
v_spgrambw(audio_in5,fs,'pJcw') ;
title('Spectrogram of Clean Speech')





