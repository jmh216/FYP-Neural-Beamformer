set(0,'DefaultAxesFontSize',13)
source_frequency = 0:10:3000;
source_angle = 0:2:360; 
figure;

load('DSB_MultiFreq_gain_6-0.05-90.mat');
subplot(332);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('6 Mic - 0.05 m Space')

load('DSB_MultiFreq_gain_6-0.15-90.mat');
subplot(335);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('6 Mic - 0.15 m Space')


load('DSB_MultiFreq_gain_6-0.25-90.mat');
subplot(338);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('6 Mic - 0.25 m Space')

% -------------------------------------------------------------------------
load('DSB_MultiFreq_gain_4-0.05-90.mat');
subplot(331);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('4 Mic - 0.05 m Space')

load('DSB_MultiFreq_gain_4-0.15-90.mat');
subplot(334);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('4 Mic - 0.15 m Space')


load('DSB_MultiFreq_gain_4-0.25-90.mat');
subplot(337);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('4 Mic - 0.25 m Space')
% -------------------------------------------------------------------------
load('DSB_MultiFreq_gain_8-0.05-90.mat');
subplot(333);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('8 Mic - 0.05 m Space')

load('DSB_MultiFreq_gain_8-0.15-90.mat');
subplot(336);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('8 Mic - 0.15 m Space')


load('DSB_MultiFreq_gain_8-0.25-90.mat');
subplot(339);surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title('8 Mic - 0.25 m Space')


% colorbar; 
% title(colorbar,'Output Gain/dB')