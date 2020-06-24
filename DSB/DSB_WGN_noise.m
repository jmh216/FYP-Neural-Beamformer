%--------------------------------------------------------------------------
% Example script for DSB response of audio input
% Author        : Jinming Hu
% Date          : 05-03-2020
% Version       : Final
%--------------------------------------------------------------------------
%% Parameter Settings
clc;
close all;
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

%% Create source and receivers

r = CreateReceiverPosition(Mic_number,Mic_spacing,Room(1),Room(2),Room(3)); % Generate receiver position [x y z] (m)
s = CreateSourcePosition(source_angle,Mic_Source_Distance,Room(1),Room(2),Room(3));

% read audio input 
filename_input = 'female_speech.wav';
[audio_in, fs] = audioread(filename_input); %sampled data, y, and a sample rate for that data, Fs=16000
t = 0:1/fs:(length(audio_in)-1)/fs;

% figure;
% plot(t,audio_in);

RIR = FindRoomImpulseResponse(c, fs, r, s, Room, beta, n);
audio_in_rir = ConvWithRIR(audio_in',RIR);

% power_input = sum(sum(abs(audio_in_rir).^2))/length(audio_in_rir)/Mic_number;
% power_input = sum(sum(abs(audio_in_rir).^2))/length(audio_in_rir)/Mic_number;
mix_in_rir = zeros(6,length(audio_in_rir));

% for i = 1:6
%     mix_in_rir(i,:) = awgn(audio_in_rir(i,:),5,'measured'); % SNR = 5 dB 
% end

%     mix_in_rir= awgn(audio_in_rir,5,'measured'); % SNR = 5 dB 

for i = 1:6
    mix_in_rir(i,:) = audio_in_rir(i,:) + rand(1,137549)/100;
end


% figure;
% plot(t, mix_in_rir);

time_length = length(mix_in_rir(1,:)); 
if mod(time_length,2) == 0 % N is even
    frequency_axis1 = (0:((time_length)/2))*Fs/(time_length); 
    frequency_axis2 = [frequency_axis1,fliplr(-frequency_axis1(2:end-1))];    
else % N is odd 
    frequency_axis1 = (0:((time_length-1)/2))*Fs/(time_length-1); 
    frequency_axis2 = [frequency_axis1,fliplr(-frequency_axis1(2:end))];
end

signal_input_freq = fft(mix_in_rir.').';  % perform fft 
weights_DSB = FindWeightsDSB(r,Mic_spacing,frequency_axis2,steering_angle,c); % calculate weights of DSB
signal_output_DSB_freq = sum(weights_DSB.*signal_input_freq); % multiply the weights and add all the channels together
signal_output_DSB_time = ifft(signal_output_DSB_freq);

figure;
subplot(221); plot(t, mix_in_rir(1,:));
xlabel("Time/s"); ylabel("Amplitude");
legend('Ref Speech Audio + WGN');
title("Waveform of noisy speech input at ref channel");

subplot(222);
v_spgrambw(mix_in_rir(1,:),fs,'pJcw') ;
title("Spectrogram of noisy speech input at ref channel");

subplot(223);
plot(t, mix_in_rir(1,:));
hold on;
plot(t,signal_output_DSB_time);
plot(t,awgn(audio_in_rir(1,:),10000000000000000000000000000000,'measured'));
hold on;
legend('Ref Noisy Audio in','DSB Audio out','Ref Target Audio in');
xlabel("Time/s"); ylabel("Amplitude");
title("Steer at 0 - Target audio at 0 - WGN");
hold off;

subplot(224);
v_spgrambw(signal_output_DSB_time,fs,'pJcw') ;
title("Spectrogram of DSB audio output");

%% Plots
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

function PlotRIR(RIR,Mic_number)
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
