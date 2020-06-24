%--------------------------------------------------------------------------
% Example script for signal genrator
% Author        : Jinming Hu
% Date          : 04-03-2020
%--------------------------------------------------------------------------
%% Parameter Settings
clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);

c = 340;                % Sound velocity (m/s)
n = 24000;               % Number of samples for RIR
beta = 0;               % Reverberation time (s)
Room = [8 8 4];			% Room dimensions [x y z] (m) 

Fs = 16000;             % Sample frequency (samples/s) 
T = 1/Fs;               % Sampling period (s)       

Mic_number = 6;         % Number of microphones
Mic_spacing = 0.15;     % Spacing of microphones (m)

source_angle = 120;
Mic_Source_Distance = 2; % Distance between source and receiver

noise_angle = 15;
Mic_Noise_Distance = 4;

%% Main

r = CreateReceiverPosition(Mic_number,Mic_spacing,Room(1),Room(2),Room(3)); % Generate receiver position [x y z] (m)
% s = [Room(1)/2+1.5 Room(2)/2+Mic_Source_Distance Room(3)/2]; % Source position [x y z] (m) 
s = CreateSourcePosition(source_angle,Mic_Source_Distance,Room(1),Room(2),Room(3));

noise_position = CreateSourcePosition(noise_angle,Mic_Noise_Distance,Room(1),Room(2),Room(3));
PlotPosition(Room,r,s,noise_position);

RIR = FindRoomImpulseResponse(c, Fs, r, s, Room, beta, n);

%% Plots
%plot rir of each channel
%--------------------------------------------------------------------------
for figure_plot = 2
figure;
for i = 1:Mic_number
    rir_plot = subplot(2,Mic_number/2,i);
    plot(RIR(i,:));
    title(rir_plot, ['Room impulse response of the ',num2str(i),'th microphone'])
    xlabel("Samples"); ylabel("Amplitude");
end
figure;
plot(RIR');
title('Room impulse response of the microphone array')
xlabel("Samples"); ylabel("Amplitude");
end
%--------------------------------------------------------------------------

%plot rir variation with change of beta
%--------------------------------------------------------------------------
% for figure_plot = 3
%     figure;
% for i = 1:4
%     rir_plot = subplot(2,2,i);
%     RIR_compare = FindRoomImpulseResponse(c, Fs, r, s, Room, i/4-0.25, n); % beta = 0-0.75
%     plot(RIR_compare(1,:));
%     title(rir_plot, ['Room impulse response of the first microphone with beta=',num2str(i/4-0.25)])
%     xlabel("Samples"); ylabel("Amplitude");
% end
% end
%--------------------------------------------------------------------------

%plot sine wave output at each channel after convolve with rir
%--------------------------------------------------------------------------
% L = 16000;              % Length of signal
% t = (0:L-1)*T;          % Time vector
% for figure_plot = 4
%     figure;
%     source_input_sine = sin(2*pi*50*t);
%     subplot(3,3,1);
%     plot(source_input_sine);
%     title('50 Hz sine wave as source')
%     xlabel("Samples"); ylabel("Amplitude");
%     
%     source_output_sine = ConvWithRIR(source_input_sine,RIR);
%     for i = 1:Mic_number
%     output_plot = subplot(3,3,i+1);
%     plot(source_output_sine(i,:));
%     title(output_plot, ['50 Hz sine wave received at the ',num2str(i),'th microphone'])
%     xlabel("Samples"); ylabel("Amplitude");
% end
% end
%--------------------------------------------------------------------------

%Read audio input and generate output file, plot signal figure
%--------------------------------------------------------------------------
for figure_plot = 5
filename_input = 'female_speech.wav';
[audio_in, fs] = audioread(filename_input); %sampled data, y, and a sample rate for that data, Fs=16000

noise_input_filename = 'ch01.wav';
[noise_in, fs_noise] = audioread(noise_input_filename);

noise_in_crop = noise_in(1:size(audio_in))*50;

t = 0:1/fs:(length(audio_in)-1)/fs;
% generate noise signal
% noise = 0.25 * sin(2*pi*1000*t);
RIR_noise =  FindRoomImpulseResponse(c, fs, r, noise_position, Room, beta, n);
audio_in_rir = ConvWithRIR(audio_in',RIR);
noise_in_rir = ConvWithRIR(noise_in_crop',RIR_noise);

figure;
t = 0:1/fs:(length(audio_in)-1)/fs;
subplot(411); plot(t,noise_in_rir); title('in(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(412); plot(t,noise_in_rir(1,:)); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(413); plot(t,noise_in_rir(6,:)); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(414); plot(t,noise_in_rir'); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
grid on;

audio_out = audio_in_rir + noise_in_rir;

figure;
t = 0:1/fs:(length(audio_in)-1)/fs;
subplot(411); plot(t,audio_out); title('in(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(412); plot(t,audio_out(1,:)); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(413); plot(t,audio_out(6,:)); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
subplot(414); plot(t,audio_out'); title('out(n)');xlabel('Time [Seconds]');ylabel('Amplitude');
grid on;

filename_output1 = 'speech_output1.wav';
audiowrite(filename_output1,audio_out(1,:),fs);
filename_output2 = 'speech_output2.wav';
audiowrite(filename_output2,audio_out(2,:),fs);
filename_output3 = 'speech_output3.wav';
audiowrite(filename_output3,audio_out(3,:),fs);
filename_output4 = 'speech_output4.wav';
audiowrite(filename_output4,audio_out(4,:),fs);
filename_output5 = 'speech_output5.wav';
audiowrite(filename_output5,audio_out(5,:),fs);
filename_output6 = 'speech_output6.wav';
audiowrite(filename_output6,audio_out(6,:),fs);
end
%--------------------------------------------------------------------------

%plot position
%--------------------------------------------------------------------------
function PlotPosition(Room,receiver_position,source_position,noise_position)
    figure;
    plot3(receiver_position(:,1),receiver_position(:,2),receiver_position(:,3),'o');
    hold on;
    plot3(source_position(:,1),source_position(:,2),source_position(:,3),'rs'); 
    plot3(noise_position(:,1),noise_position(:,2),noise_position(:,3),'x');
    legend('Microphones','Signal Source','Noise');
    axis([0 Room(1) 0 Room(2) 0 Room(3)]);
    xlabel("x-axis/m"); ylabel("y-axis/m"); zlabel("z-axis/m");
    title("Position of source and microphone array");
    grid on; box on; axis square;hold off;
end
%--------------------------------------------------------------------------

%% Functions
function receiver_position = CreateReceiverPosition(Mic_number,Micro_Spacing,x,y,z)
    receiver_position = zeros(Mic_number,3); % A matrix with size nMic x 3, 3 for three axis x ,y, z
    Array_length = Micro_Spacing*(Mic_number-1);
    for index = 1:Mic_number
        receiver_position(index,:) = [(x-Array_length)/2+(index-1)*Micro_Spacing y/2 z/2]; % Microphones'position different in x axis
    end
end

function source_position = CreateSourcePosition(source_degree,radius,x,y,z) 
    source_position = zeros(1,3);
    source_degree_radian = source_degree*pi/180;
    source_position(1,:) = [x/2+radius*cos(source_degree_radian) y/2+radius*sin(source_degree_radian) z/2]; 
end

function rir_array = FindRoomImpulseResponse(c, fs, receiver_position, source_position, Room, beta, n) 
    Mic_number = size(receiver_position,1);    
    rir_array = zeros(Mic_number,n);
    for index = 1:Mic_number
        h = rir_generator(c, fs, receiver_position(index,:), source_position, Room, beta, n); 
        rir_array(index,:) = h;
    end
end

function signal_output = ConvWithRIR(signal_input,RIR) 
    Mic_number = size(RIR,1); 
    sample_number = size(signal_input,2);
    signal_output = zeros(Mic_number,sample_number);
    for index = 1:Mic_number
        signal_output(index,:) = filter(RIR(index,:),1,signal_input); 
    end
end


