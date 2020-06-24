%% setup paths
restoredefaultpath
addpath('submodules/sap-elobes-microphone-arrays');
addpath('submodules/sap-elobes-utilities');
addpath('submodules/sap-voicebox/voicebox');
addpath('lib');

clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);

%% beamformer design
% ensure empty struct
spec = struct();

% sample rate
spec.fs = 16000;            

% cartesian unit vector
spec.look_unit_vec = [1 0 0];  % default
% spec.look_unit_vec = [0 1 0];  

% microphone array
nSensors = 6;
spacing = 0.15;
spec.ema_fcn = @()FreeField_20180202_01_ULA_configurable(nSensors,spacing);        % function handle for elobes microphone array to use 
                                                                                   % create a sensor position matrix 
filtCoefs = design_beamformer(spec);

figure;
plot(filtCoefs)

%% demonstrate effect on signal from look direction
ema = spec.ema_fcn();
ema.setSoundSpeed(v_soundspeed);
ema.prepareData(spec.fs);
[target_az,target_inc,~] = mycart2sph(spec.look_unit_vec);
h = ema.getImpulseResponseForSrc(target_az,target_inc);

figure;
subplot(411);plot(h);title('Array response')
subplot(412);plot(filtCoefs);title('Beamformer coefficients')
subplot(413);plot(fftfilt(filtCoefs,[h;zeros(size(filtCoefs,1)-1,nSensors)]));title('Per channel filtered response')
subplot(414);plot(sum(fftfilt(filtCoefs,[h;zeros(size(filtCoefs,1)-1,nSensors)]),2),'linewidth',1.5);title('Beamformer output')
linkaxes(get(gcf,'children'))

%% example audio
[speech,fs] = v_readwav('./submodules/sap-voicebox/voicebox_demo/data/as01b0.wav');
% [speech,fs] = v_readwav('female_speech.wav');
mic_speech = fftfilt(h,speech);

diffuse_noise = zeros(size(mic_speech));
noise_az_deg = (0:5:355).';
h_circ = ema.getImpulseResponseForSrc(deg2rad(noise_az_deg),deg2rad(90*ones(size(noise_az_deg))));
for idoa = 1:length(noise_az_deg)
    new_mono_noise = v_stdspectrum(11,'t',fs,length(speech));
    diffuse_noise = diffuse_noise + 1/10000 * fftfilt(squeeze(h_circ(:,:,idoa)),new_mono_noise);
end

mixed = mic_speech + diffuse_noise;

mixed_ref = mixed(:,ema.refChan);
target_out = sum(fftfilt(filtCoefs,mic_speech),2);
bf_out = sum(fftfilt(filtCoefs,mixed),2);

figure;
plot(mixed_ref);hold all;plot(bf_out);plot(target_out);
legend('mixed_{ref}','bf_{out}','target_{out}')
    
