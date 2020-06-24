clc;
close all;
clear variables;
set(0,'DefaultAxesFontSize',11)
set(0,'DefaultLineLineWidth',1);
%% setup paths
% restoredefaultpath
% restoredefaultpath
addpath('submodules/sap-elobes-microphone-arrays');
addpath('submodules/sap-elobes-utilities');
addpath('submodules/sap-voicebox/voicebox');
addpath('lib');

Audio_in_folder = 'C:\Users\74778\Desktop\Audio files\et05_ped_real';
Audio_out_folder_DSB = 'C:\Users\74778\Desktop\Audio files DSB\et05_ped_real';
Audio_out_folder_MVDR = 'C:\Users\74778\Desktop\Audio files MVDR\et05_ped_real';
% Audio_out_folder_clean = 'C:\Users\74778\Desktop\Audio files clean\et05_ped_real';
% Audio_out_folder_ref = 'C:\Users\74778\Desktop\Audio files ref\et05_ped_real';

if ~isfolder(Audio_in_folder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', Audio_in_folder);
    uiwait(warndlg(errorMessage));
    Audio_in_folder = uigetdir(); % Ask for a new one.
    if Audio_in_folder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(Audio_in_folder, '*.wav'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for j = 1 : length(theFiles)/7
    channel = [];
    bf_out = [];
    for k = 7*j-6 : 7*j
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        fprintf(1, 'Now reading %s\n', fullFileName);
        [channel(:,k-(j-1)*7), fs] =  v_readwav(fullFileName);    
    end
    multi_input = channel(:,[2:7]);
    ch0_filename = theFiles(7*j-6).name;
    ch5_filename = theFiles(7*j-6+5).name;

    %% beamformer design
    % ensure empty struct
    spec = struct();
    % sample rate
    spec.fs = 16000;            
    % Create cartesian unit vector for looking direction
    steer_angle = 90; % 90 degree-[0 1 0] ; 0 degree-[1 0 0]
    unit_source_position = [cosd(steer_angle) sind(steer_angle) 0]; 
    spec.look_unit_vec = unit_source_position;  
    % microphone array
    nSensors = 6;
    spacing = 0.05;
    spec.ema_fcn = @()FreeField_20180202_01_ULA_configurable(nSensors,spacing);        % function handle for elobes microphone array to use % create a sensor position matrix                                                                    

    %% example audio
    t = 0:1/fs:(length(multi_input)-1)/fs;

    for i = 1:2
        if i ==1
            spec.noise_model = 'SPH_ISO';
            filtCoefs = design_beamformer(spec);
        end
        if i == 2
            spec.noise_model = 'identity';
            filtCoefs = design_beamformer(spec);
        end
        ema = spec.ema_fcn();
        ema.setSoundSpeed(v_soundspeed);
        ema.prepareData(spec.fs);

        mic_speech = multi_input;

        mixed_ref = mic_speech(:,ema.refChan);
        bf_out(:,i) = sum(fftfilt(filtCoefs,mic_speech),2);
    end

%     figure;
%     subplot(121);
%     plot(t,mixed_ref);hold all;plot(t,bf_out(:,1));plot(t,channel(:,1));
%     xlabel("Time/s"); ylabel("Amplitude");
%     legend('mixed reference','MVDR BF Directional out','target out')
%     title('MVDR audio result (Directional noise model）');
%     hold off;
%     
%     subplot(122);
%     plot(t,mixed_ref);hold on;plot(t,bf_out(:,2));plot(t,channel(:,1));
%     xlabel("Time/s"); ylabel("Amplitude");
%     legend('mixed reference','DSB out','target out');
%     title('DSB audio result');
%     
%     figure;
%     subplot(121);
%     v_spgrambw(bf_out(:,1),fs,'pJcw') ;
%     title('Spectrogram of MVDR audio result (Directional noise model）');
%     subplot(122);
%     v_spgrambw(bf_out(:,2),fs,'pJcw') ;
%     title('Spectrogram of DSB audio result');

    underlineLocations = find(ch0_filename == '.');
    curr1 = ch0_filename(1:underlineLocations(1)-1);
    filename = curr1;

    % [folder, baseFileName, extension] = fileparts(filename);
    
%     % ref channel    
%     underlineLocations = find(ch5_filename == '.');
%     curr5 = ch5_filename(1:underlineLocations(1)-1);
%     
%     filename_ref = fullfile(Audio_out_folder_ref, [curr5,'_ref.wav']);
%     audiowrite(filename_ref,multi_input(:,5),fs);
    
%     % clean
%     filename_clean = fullfile(Audio_out_folder_clean, [filename,'_clean.wav']);
%     audiowrite(filename_clean,channel(:,1),fs);
% 
    filename_MVDR = fullfile(Audio_out_folder_MVDR, [filename,'_MVDR.wav']);
    audiowrite(filename_MVDR,bf_out(:,1),fs);

    filename_DSB = fullfile(Audio_out_folder_DSB, [filename,'_DSB.wav']);
    audiowrite(filename_DSB,bf_out(:,2),fs);
    
end

%% Functions
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
        noise_angle = 120; % 90 degree-[0 1 0] ; 0 degree-[1 0 0]
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
         R = zeros(size(d_t,1),size(d_t,1),size(d_t,2));
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

