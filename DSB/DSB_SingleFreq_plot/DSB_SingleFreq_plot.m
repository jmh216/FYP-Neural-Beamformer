load('DSB_Single_gain_6mic-0.15spacing-90-500Hz.mat');
signal_output_DSB_gain_500 = signal_output_DSB_gain;
load('DSB_Single_gain_6mic-0.15spacing-90-1000Hz.mat');
signal_output_DSB_gain_1000 = signal_output_DSB_gain;
load('DSB_Single_gain_6mic-0.15spacing-90-2000Hz.mat');
signal_output_DSB_gain_2000 = signal_output_DSB_gain;
load('DSB_Single_gain_6mic-0.15spacing-90-3000Hz.mat');
signal_output_DSB_gain_3000 = signal_output_DSB_gain;

source_angle = 0:1:360; 


figure;
subplot(221);plot(source_angle,10*log10(signal_output_DSB_gain_500));xlabel("Angle/degree"); ylabel("Gain/dB");
title('DSB Output Gain of 500 Hz sinusoidal signal'),xlim([0 360]);
subplot(222);plot(source_angle,10*log10(signal_output_DSB_gain_1000));xlabel("Angle/degree"); ylabel("Gain/dB");
title('DSB Output Gain of 1000 Hz sinusoidal signal'),xlim([0 360]);
subplot(223);plot(source_angle,10*log10(signal_output_DSB_gain_2000));xlabel("Angle/degree"); ylabel("Gain/dB");
title('DSB Output Gain of 2000 Hz sinusoidal signal'),xlim([0 360]);
subplot(224);plot(source_angle,10*log10(signal_output_DSB_gain_3000));xlabel("Angle/degree"); ylabel("Gain/dB");
title('DSB Output Gain of 3000 Hz sinusoidal signal'),xlim([0 360]);

figure;
subplot(221);
output_gain_dB = convert_dB(signal_output_DSB_gain_500);
h1 = polarplot(source_angle*pi/180,output_gain_dB+40,'Linewidth',1);
haxes1 = get(h1,'Parent');
haxes1.RTickLabel = {'-40','-30','-20','-10','0 dB'};
title('Polar plot of output gain of 500 Hz sinusoidal signal')

subplot(222);
output_gain_dB = convert_dB(signal_output_DSB_gain_1000);
h2 = polarplot(source_angle*pi/180,output_gain_dB+40,'Linewidth',1);
haxes2 = get(h2,'Parent');
haxes2.RTickLabel = {'-40','-30','-20','-10','0 dB'};
title('Polar plot of output gain of 1000 Hz sinusoidal signal')

subplot(223);
output_gain_dB = convert_dB(signal_output_DSB_gain_2000);
h3 = polarplot(source_angle*pi/180,output_gain_dB+40,'Linewidth',1);
haxes3 = get(h3,'Parent');
haxes3.RTickLabel = {'-40','-30','-20','-10','0 dB'};
title('Polar plot of output gain of 2000 Hz sinusoidal signal')

subplot(224);
output_gain_dB = convert_dB(signal_output_DSB_gain_3000);
h4 = polarplot(source_angle*pi/180,output_gain_dB+40,'Linewidth',1);
haxes4 = get(h4,'Parent');
haxes4.RTickLabel = {'-40','-30','-20','-10','0 dB'};
title('Polar plot of output gain of 3000 Hz sinusoidal signal')

function output_gain_dB = convert_dB(signal_output_DSB_gain)
output_gain_dB = 10*log10(signal_output_DSB_gain);
output_gain_dB = output_gain_dB - max(output_gain_dB);
output_gain_dB(output_gain_dB<-40) = -40;
end