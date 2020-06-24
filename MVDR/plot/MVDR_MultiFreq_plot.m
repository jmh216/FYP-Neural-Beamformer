close all;
clc;

set(0,'DefaultAxesFontSize',13)
source_frequency = 0:30:3000;
source_angle = 0:5:360; 

% ----------------------------------------------------------------------- %
figure;
load('4-0.05-0-iso.mat');
subplot(221);
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 4 Mic 0.05 m Steer 0','Spherical isotropic noise'})

load('6-0.05-0-iso.mat');
subplot(222);
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 6 Mic 0.05 m Steer 0','Spherical isotropic noise'})


load('4-0.15-0-iso.mat');
subplot(223);
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 4 Mic 0.15 m Steer 0','Spherical isotropic noise'})

load('6-0.15-0-iso.mat');
subplot(224);
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 6 Mic 0.15 m Steer 0','Spherical isotropic noise'})

% ----------------------------------------------------------------------- %
figure;
subplot(121);
load('6-0.05-0-identity.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 6 Mic 0.05 m Steer 0','White noise'})

subplot(122);
load('6-0.15-0-identity.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz') 
title({'MVDR 6 Mic 0.15 m Steer 0','White noise'})
% ----------------------------------------------------------------------- %

figure;
subplot(221);
load('6-0.05-0-45-dirc.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz');
title({'MVDR 6 Mic 0.05 m Steer 0','Directional noise from 45'})

subplot(222);
load('6-0.05-0-90-dirc.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz');
title({'MVDR 6 Mic 0.05 m Steer 0','Directional noise from 90'})

subplot(223);
load('6-0.05-0-150-dirc.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz');
title({'MVDR 6 Mic 0.05 m Steer 0','Directional noise from 150'})

subplot(224);
load('6-0.05-0-180-dirc.mat');
surf(source_angle,source_frequency,10*log10(final_out_gain))
view(0,90) 
xlim([0 360]);
shading interp	%remove grid
xlabel('Angle/degree'); 
ylabel('Frequency/Hz');
title({'MVDR 6 Mic 0.05 m Steer 0','Directional noise from 180'})


% colorbar; 
% title(colorbar,'Output Gain/dB')