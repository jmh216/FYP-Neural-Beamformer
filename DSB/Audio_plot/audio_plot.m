

filename_input = 'female_speech.wav';
[audio_in, fs] = audioread(filename_input); %sampled data, y, and a sample rate for that data, Fs=16000
t = 0:1/fs:(length(audio_in)-1)/fs;

load('6-0.15-2-0-0.mat');
signal_output_0 = signal_output_DSB_time;
load('6-0.15-2-0-45.mat');
signal_output_45 = signal_output_DSB_time;
load('6-0.15-2-0-90.mat');
signal_output_90 = signal_output_DSB_time;
load('6-0.15-2-0-180.mat');
signal_output_180 = signal_output_DSB_time;

figure;
subplot(221);
plot(t,audio_in_rir(1,:));
hold on;
plot(t,signal_output_0);
xlabel("Time/s"); ylabel("Amplitude");
legend('Ref Audio in','DSB Audio out')
title('Steer at 0 - Audio input at 0 degree');

subplot(222);
plot(t,audio_in_rir(1,:));
hold on;
plot(t,signal_output_45);
xlabel("Time/s"); ylabel("Amplitude");
legend('Ref Audio in','DSB Audio out')
title('Steer at 0 - Audio input at 45 degree');

subplot(223);
plot(t,audio_in_rir(1,:));
hold on;
plot(t,signal_output_90);
xlabel("Time/s"); ylabel("Amplitude");
legend('Ref Audio in','DSB Audio out')
title('Steer at 0 - Audio input at 90 degree');

subplot(224);
plot(t,audio_in_rir(1,:));
hold on;
plot(t,signal_output_180);
xlabel("Time/s"); ylabel("Amplitude");
legend('Ref Audio in','DSB Audio out')
title('Steer at 0 - Audio input at 180 degree');

% load('6-0.15-2-90-0.mat');
% signal_output_0 = signal_output_DSB_time;
% load('6-0.15-2-90-30.mat');
% signal_output_30 = signal_output_DSB_time;
% load('6-0.15-2-90-90.mat');
% signal_output_90 = signal_output_DSB_time;
% load('6-0.15-2-90-120.mat');
% signal_output_120 = signal_output_DSB_time;
% 
% 
% figure;
% subplot(221);
% plot(t,audio_in/25);
% hold on;
% plot(t,signal_output_0);
% xlabel("Time/s"); ylabel("Amplitude");
% legend('Audio input','Audio output')
% title('Steer at 90 with input audio at 0 degree');
% 
% subplot(222);
% plot(t,audio_in/25);
% hold on;
% plot(t,signal_output_30);
% xlabel("Time/s"); ylabel("Amplitude");
% legend('Audio input','Audio output')
% title('Steer at 90 with input audio at 30 degree');
% 
% subplot(223);
% plot(t,audio_in/25);
% hold on;
% plot(t,signal_output_90);
% xlabel("Time/s"); ylabel("Amplitude");
% legend('Audio input','Audio output')
% title('Steer at 90 with input audio at 90 degree');
% 
% subplot(224);
% plot(t,audio_in/25);
% hold on;
% plot(t,signal_output_120);
% xlabel("Time/s"); ylabel("Amplitude");
% legend('Audio input','Audio output')
% title('Steer at 90 with input audio at 120 degree');