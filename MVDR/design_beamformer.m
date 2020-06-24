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
end

% design the filter in frequency domain
% (include the conjugate so we can apply in time domain using convolution)
wl_tmp = conj(mvdr_fcn(R, d_t, params.mvdr_debug_level).'); %[nFreq,nChans]

% get fir coefficients
% (do real inverse fft and shift to make it causal)
coefs = circshift(irfft(wl_tmp,nfft,1),floor(nfft/2),1);


