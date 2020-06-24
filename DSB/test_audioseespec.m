clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);

filename_output1 = '1.wav';
    [audio_in1, fs] = audioread(filename_output1); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output2 = '2.wav';
    [audio_in2, fs] = audioread(filename_output2); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output3 = '3.wav';
    [audio_in3, fs] = audioread(filename_output3); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output4 = '4.wav';
    [audio_in4, fs] = audioread(filename_output4); %sampled data, y, and a sample rate for that data, Fs=16000
filename_output5 = 'ori.wav';
    [audio_in5, fs] = audioread(filename_output5); %sampled data, y, and a sample rate for that data, Fs=16000
    
t = 0:1/fs:(length(audio_in1)-1)/fs;

figure;
subplot(221);
plot(t,audio_in1);
hold on;
plot(t,audio_in5)
legend('REF','Clean');
hold off;

subplot(222);
plot(t,audio_in2);
hold on;
plot(t,audio_in5)
legend('DSB','Clean');
hold off;

subplot(224);
plot(t,audio_in3([1:length(t)],1));
hold on;
plot(t,audio_in5)
legend('BLSTM+GEV+BAN','Clean');
hold off;

subplot(223);
plot(t,audio_in4);
hold on;
plot(t,audio_in5)
legend('MVDR','Clean');
title('Waveform Comparison')
hold off;

figure;
title('Spectrogram Comparison');

subplot(221);
v_spgrambw(audio_in1,fs,'pJcw') ;
title('REF')

subplot(222);
v_spgrambw(audio_in2,fs,'pJcw') ;
title('DSB')


subplot(224);
v_spgrambw(audio_in3,fs,'pJcw') ;
title('BLSTM+GEV+BAN')

subplot(223);
v_spgrambw(audio_in4,fs,'pJcw') ;
title('MVDR')

figure;
v_spgrambw(audio_in5,fs,'pJcw') ;
title('Spectrogram of Clean Speech')





