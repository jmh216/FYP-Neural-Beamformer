clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);
%% setup paths
% restoredefaultpath
addpath('submodules/sap-elobes-microphone-arrays');
addpath('submodules/sap-elobes-utilities');
addpath('submodules/sap-voicebox/voicebox');
addpath('lib');

%% beamformer design
% ensure empty struct
spec = struct();

% sample rate
spec.fs = 16000;            

% cartesian unit vector
% spec.look_unit_vec = [1 0 0];  % oringinal 180 degree
% spec.look_unit_vec = [0 1 0];  

% Create cartesian unit vector for looking direction
steer_angle = 90; % 90 degree-[0 1 0] ; 0 degree-[1 0 0]
unit_source_position = [cosd(steer_angle) sind(steer_angle) 0]; 
spec.look_unit_vec = unit_source_position;  

% microphone array
nSensors = 6; % not useful when testing, fixed config
spacing = 0.05; % not useful when testing, fixed config
spec.ema_fcn = @()FreeField_20180202_01_ULA_configurable(nSensors,spacing);        % function handle for elobes microphone array to use % create a sensor position matrix                                                                    
spec.noise_model = 'SPH_ISO';
% spec.noise_model = 'directional';
% spec.noise_model = 'identity';
% spec.noise_model = 'multi-directional';

filtCoefs = design_beamformer(spec);

% figure;
% plot(filtCoefs)

%% Beamformer response

source_degree = 360;    % Angle of source need to plot
source_resolution = 1;  % Angle difference of each source
source_radius = 1;      % Radius of the source to the center of the array Default:2
s = CreateSourcePosition(source_degree,source_resolution, source_radius, 0, 0, 0); % Generate source position [x y z] (m)
source_angle = 0:source_resolution:source_degree;  

% source_angle = 0:source_resolution:source_degree; 
% source_input = cos(2*pi*source_frequency*t); 

output_gain = zeros(1,size(s,1));

for source_index = 1:size(s,1) % for each source
    ema = spec.ema_fcn();
    ema.setSoundSpeed(v_soundspeed);
    ema.prepareData(spec.fs);
    
    Fs = 16000;             % Sample frequency (samples/s) 
    T = 1/Fs;               % Sampling period (s)       
    L = 1600;               % Length of signal
    t = (0:L)*T;            % Time vector
    
    % sinusoidal input
    source_frequency = 1000;
    source = cos(2*pi*source_frequency*t); 
    
    % generate RIR for this source     
    [target_az,target_inc,~] = mycart2sph(s(source_index,:));
    h = ema.getImpulseResponseForSrc(target_az,target_inc);
    
%   fft
    source_input = fftfilt(h,source);
%   ------------------------- add isotrophic noise ------------------------
    diffuse_noise = zeros(size(source_input));
    noise_az_deg = (0:5:355).';
    h_circ = ema.getImpulseResponseForSrc(deg2rad(noise_az_deg),deg2rad(90*ones(size(noise_az_deg))));
    for idoa = 1:length(noise_az_deg)
        new_mono_noise = v_stdspectrum(11,'t',Fs,length(source_input));
        diffuse_noise = diffuse_noise + 1/10000 * fftfilt(squeeze(h_circ(:,:,idoa)),new_mono_noise);
    end
% -------------------------------------------------------------------------
    % add diffuse noise
%     mixed = source_input + diffuse_noise;
    
    % do not add diffuse noise
    mixed = source_input;
    
    %MVDR BF
    bf_out = sum(fftfilt(filtCoefs,mixed),2);
    
    
%     power_input = sum(sum(abs(mixed).^2))/length(source_input)/nSensors;
    
    power_input = sum(sum(abs(source_input).^2))/length(source_input)/nSensors;
    power_output = sum(sum(abs(bf_out).^2))/length(bf_out);
    output_gain(source_index) = (power_output/power_input); 
% -------------------------------------------------------------------------
end

figure;
plot3(ema.sensorCartesianPositionsDefault(:,1),ema.sensorCartesianPositionsDefault(:,2),ema.sensorCartesianPositionsDefault(:,3),'o');
hold on;
plot3(unit_source_position(1),unit_source_position(2),unit_source_position(3),'x');
axis([-2 2 -2 2 -2 2]);
xlabel("x-axis/m"); ylabel("y-axis/m"); zlabel("z-axis/m");
legend('Microphone','Source')
title("Position of microphone array");
grid on; box on; axis square;hold off;

figure;
subplot(121);
plot(source_angle,10*log10(output_gain));
xlabel("Angle/degree"); ylabel("Gain/dB");
% title({["MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Directional noise 45 degree" ]});
% title({["MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Spherical isotropic noise"]});
% title({["MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - White noise"]});
title({["For array geometry used in testing:"],[" MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Spherical isotropic noise"]});
% title({["For array geometry used in testing:"],[" DSB beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree"]});

xlim([0 360]);
subplot(122);
output_gain_dB = 10*log10(output_gain);
output_gain_dB = output_gain_dB - max(output_gain_dB);
output_gain_dB(output_gain_dB<-40) = -40;
h = polarplot(source_angle*pi/180,output_gain_dB+40,'Linewidth',1);
haxes = get(h,'Parent');
haxes.RTickLabel = {'-40','-30','-20','-10','0 dB'};
% title({["Polar plot of MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Directional noise 45 degree"]});
% title({["Polar plot of MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Spherical isotropic noise"]});
% title({["Polar plot of MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - White noise"]});
title({["For array geometry used in testing:"],[" Polar plot of MVDR beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree - Spherical isotropic noise"]});
% title({["For array geometry used in testing:"],[" Polar plot of DSB beampattern for 1000 Hz signal" ],["Look direction " + num2str(steer_angle) + " degree"]});

%% Functions
function source_position = CreateSourcePosition(source_degree,source_resolution,radius,x,y,z) 
    source_position = zeros(source_degree/source_resolution+1,3); % eg. 360/60+1, 1 more for 0 degree
    source_angle = 0:source_resolution*pi/180:source_degree*pi/180;
    for index = 1:source_degree/source_resolution+1
        source_position(index,:) = [x/2+radius*cos(source_angle(index)) y/2+radius*sin(source_angle(index)) z/2]; 
    end
end

function[coefs] = design_beamformer(in)
% calculates filter coefficients according to the design specifications
% given as fields of a struct
%

% Required parameters (with no default assumption) specified here as empty fields
defaults.fs = [];             % sample rate
defaults.look_unit_vec = [];  % cartesian unit vector
defaults.ema_fcn = [];        % function handle for elobes microphone array to use 

% Optional parameters which are assumned unless overridden

defaults.noise_model = 'SPH_ISO'; % TODO: Implement more options (e.g. 'WHITE','FROM_FILE')
% defaults.noise_model = 'directional'; % TODO: Implement more options (e.g. 'WHITE','FROM_FILE')
% defaults.noise_model = 'identity'; % TODO: Implement more options (e.g. 'WHITE','FROM_FILE')

defaults.soundspeed = v_soundspeed();
defaults.refChanName = 'refChan'; % property of the microphone array
defaults.mvdr_debug_level = 0;      % can opt to cause warnings or errors on singularity
defaults.max_condition_number = 100;


% Optional parameters which are not assumed
defaults.noise_example_path = ''; % TODO: Allow a file to be specified
defaults.refChanNumber = [];
defaults.fir_len = [];            %


% Hard coded parameters
mvdr_fcn = @fcn_20170222_01_mvdr_weights; % function used to calculate weights

% validate mutually exclusive options
if isfield(in,'refChanName') & isfield(in,'refChanNumber')
    error('Only specify one of refChanName and refChanNumber')
end

if isfield(in,'noise_example_path') & ~isequal(in.noise_model,'FROM_FILE')
    warning('noise_example_path is specified but unused. Set noise_model=FROM_FILE to use it.')
end

% 
params = override_valid_fields(defaults,in);


%% validation

% sample rate
if ~isscalar(params.fs) || ~isreal(params.fs)
    error('fs must be a real scalar')
end

% look direction
sz_look = size(params.look_unit_vec);
if sz_look ~= [1,3]
    if sz_look == [3,1]
        % ok we can deal with that
        params.look_unit_vec = params.look_unit_vec.';
    else
        error('look_unit_vec should be a [1,3] vector of cartesian coordinates')
    end
end
if cartnorm(params.look_unit_vec)~=1
    error('look_unit_vec should be a unit vector')
end

% microphone array - not much to be done... just try to instatiate it
ema = params.ema_fcn();

% get channel number for reference
if ~isempty(params.refChanNumber)
    refChan = params.refChanNumber;
else
    refChan = ema.(params.refChanName);
end

ema.setSoundSpeed(params.soundspeed);
ema.prepareData(params.fs);


% get steering vector
[az,inc,r] = mycart2sph(params.look_unit_vec);
h = ema.getImpulseResponseForSrc(az,inc);
ir_crop = size(h,1); % start with uncropped response
nfft = ir_crop;

% adjust if required for specified fir_len
if ~isempty(params.fir_len)
    ir_crop = min(ir_crop,params.fir_len);
    nfft = max(nfft,params.fir_len);
end

% get RTF
H = rfft(h(1:ir_crop,:),nfft,1);
d = H./H(:,refChan); % rtf with respect to reference channel
d_t = d.'; %[nCham, nFreq]

% get noise covariance
switch params.noise_model
    case 'SPH_ISO'
        % ideally doas would be equally spaced on sphere. Approximate it...
        load('sph_dist_794.mat','az_794','inc_794');
        h = ema.getImpulseResponseForSrc(az_794,inc_794);
        H = rfft(h(1:ir_crop,:,:),nfft,1);
        [nFreq,nChan,nDOA] = size(H);

        % covariance for each direction
        R = bsxfun(@times,permute(H,[2 4 1 3]),conj(permute(H,[4 2 1 3]))); %[nChan nChan nFreq nDOA]
        % covariance assuming same power from each direction
        R = mean(R,4); %[nChan nChan nFreq]
        
        % at low frequencies especially Riso is rank deficient so regularise
        for ifreq = 1:nFreq
            R(:,:,ifreq) = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                R(:,:,ifreq),params.max_condition_number);
        end
    
    case 'directional'           
        noise_angle = 45; % 90 degree-[0 1 0] ; 0 degree-[1 0 0]
        unit_noise_position = [cosd(noise_angle) sind(noise_angle) 0]; 
        noise_cart = unit_noise_position;  
        [noise_az,noise_inc,~] = mycart2sph(noise_cart);
        
        h = ema.getImpulseResponseForSrc(noise_az,noise_inc);
%         h = ema.getImpulseResponseForSrc(deg2rad(70),deg2rad(0));

        H = rfft(h(1:ir_crop,:,:),nfft,1);
        [nFreq,nChan,nDOA] = size(H);

        % covariance for each direction
        R = bsxfun(@times,permute(H,[2 4 1 3]),conj(permute(H,[4 2 1 3]))); %[nChan nChan nFreq nDOA]

        % at low frequencies especially Riso is rank deficient so regularise
        for ifreq = 1:nFreq
            R(:,:,ifreq) = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                R(:,:,ifreq),params.max_condition_number);
        end
        
     case 'identity'    
         R=zeros(size(d_t,1),size(d_t,1),size(d_t,2));
         for i = 1:size(d_t,2)
            R(:,:,i) = eye(size(d_t,1));
         end
        
    case 'multi-directional'    
            freq_range = 10;
            start_frrq = 90;
            noise_az_multi = zeros(freq_range+1);
            noise_inc_multi = zeros(freq_range+1);
            for noise_angle = start_frrq : start_frrq+freq_range
                unit_noise_position = [cosd(noise_angle) sind(noise_angle) 0]; 
                noise_cart = unit_noise_position;  
                [noise_az_multi(noise_angle-start_frrq+1),noise_inc_multi(noise_angle-start_frrq+1),~] = mycart2sph(noise_cart);             
            end
            
            h = ema.getImpulseResponseForSrc(noise_az_multi,noise_inc_multi);
            H = rfft(h(1:ir_crop,:,:),nfft,1);
            [nFreq,nChan,nDOA] = size(H);

            % covariance for each direction
            R = bsxfun(@times,permute(H,[2 4 1 3]),conj(permute(H,[4 2 1 3]))); %[nChan nChan nFreq nDOA]

            % at low frequencies especially Riso is rank deficient so regularise
            for ifreq = 1:nFreq
                R(:,:,ifreq) = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                    R(:,:,ifreq),params.max_condition_number);
            end
end

% design the filter in frequency domain
% (include the conjugate so we can apply in time domain using convolution)
wl_tmp = conj(mvdr_fcn(R, d_t, params.mvdr_debug_level).'); %[nFreq,nChans]


% get fir coefficients
% (do real inverse fft and shift to make it causal)
coefs = circshift(irfft(wl_tmp,nfft,1),floor(nfft/2),1);
end

