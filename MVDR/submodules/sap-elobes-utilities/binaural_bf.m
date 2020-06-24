function y = binaural_bf(sig,air,pp)
% BLIND_BF performs blind beamforming using multichannel noisy inputs. Both the
% relative transfer function (RTF) and the noise covariance matrix are
% estimated blindly based on the existing methods.
%
% INPUTS:
% sig:  time-domain multichannel signal, with one column as one channel.
%     
% pp:   parameters of the algorithm. Default values are shown in brackets.
%       pp.win_len(256): Window size of short-time-Fourier-transform
%       (STFT). 
%       pp.overlap(0.5): Overlap between frames in range [0,1]
%       pp.f_debug(0):   Whether operates in the debug mode 
%       
%       pp.beamformer('mwf'): choose the beamformer type
%       pp.noi_est('s'): If 's', use souden's method [4] for noise
%       covariance matrix estimation. If 'h', use [3] instead.
%       pp.channel(1): the reference channel
% OUTPUTS:
% y:    time-domain single-channel beamformer output
% w_mat:n_mic*n_fft_2*n_frame matrix consisting of frequency-domain filter
%       coefficients in each frame.
%
% REFERENCES: 
% [1] Benesty, Jacob, Jingdong Chen, and Yiteng Huang.
% Microphone Array Signal Processing. Springer Science & Business Media,
% 2008.
% [2] Pan, Chao, Jingdong Chen, and J. Benesty. Performance Study
% of the MVDR Beamformer as a Function of the Source Incidence Angle.
% IEEE/ACM Transactions on Audio, Speech, and Language Processing 22, no. 1
% (2014): 67-79.
% [3] Hendriks, R.C., and T. Gerkmann. Noise Correlation Matrix Estimation for
% Multi-Microphone Speech Enhancement. IEEE Transactions on Audio, Speech, and Language 
% Processing 20, no. 1 (2012): 223-233. doi:10.1109/TASL.2011.2159711.
% [4] Souden, M., Jingdong Chen, J. Benesty, and S. Affes. an Integrated Solution for Online Multichannel Noise Tracking and Reduction.
% IEEE Transactions on Audio, Speech, and Language Processing 19, no. 7 (September 2011).
% 
% Wei Xue, Imperial College, Feb 14, 2017


%--- Check and unpack parameters ---
if nargin<2
    pp = [];
end

if isstruct(pp) || isempty(pp)
    default_fields = {'win_len','overlap','f_debug','f_rtf','smooth_factor'};
    default_values = [256, 0.5, 0, 0,0.8];
    n_fields = length(default_fields);
	for i_field = 1:n_fields
		if ~isfield(pp,default_fields{i_field})
			eval([default_fields{i_field} '=' num2str(default_values(i_field)) ';']);
		else
			eval([default_fields{i_field} '=pp.' default_fields{i_field} ';']);
		end
	end

    if ~isfield(pp,'beamformer')
        pp.beamformer = 'mvdr';
    end

else
	error('please specify the parameters.');	
end

sig_f_multi = rfft_multi(sig,win_len,overlap);
[n_mic,n_fft_2,n_frame] = size(sig_f_multi);


%--------------------perform mvdr/wmf------------------------
%Initialize ryy and rvv, set parameters
ryy_mat = stft2rxx(sig_f_multi,smooth_factor);
% n_init = fix(pp.fs/(win_len*(1-overlap))/2);
rvv = mean(ryy_mat(:,:,:,1:50),4);
spp = zeros(n_fft_2,1);
qq = struct('spp',spp,'phi_0',-2e3);    

rtfout = rir2stv(air.h,n_fft_2);

clean_sig_f = zeros(2,n_fft_2,n_frame);


minVal = 10^(-25/10);
spp_mat = zeros(n_fft_2,n_frame);

%framewise processing
for i_frame = 1:n_frame    
    sig_f = sig_f_multi(:,:,i_frame);
    ryy = ryy_mat(:,:,:,i_frame);
    qq.ryy = ryy;
    [rvv,qq] = estnoicov_s(sig_f,rvv,qq);
    spp_mat(:,i_frame) = qq.spp;

    
    %_______________beamforming_______________
    for i_freq = 1:n_fft_2
        h_bf = h_mvdr(rvv(:,:,i_freq),rtfout(:,i_freq)./rtfout(1,i_freq)); %mvdr
        clean_sig_f(1,i_freq,i_frame) = h_bf'*sig_f(:,i_freq); %spatial filtering        
        
        if strcmp(pp.beamformer,'mwf')
            pyyo = h_bf'* ryy(:,:,i_freq)*h_bf;
            pvvo = h_bf'* rvv(:,:,i_freq)*h_bf;
            clean_sig_f(1,i_freq,i_frame) = clean_sig_f(1,i_freq,i_frame).*max(minVal,1-pvvo/pyyo);            
        end

        clean_sig_f(2,i_freq,i_frame) = clean_sig_f(1,i_freq,i_frame)*rtfout(2,i_freq)./rtfout(1,i_freq);
    end   
end

y = irfft_multi(clean_sig_f,win_len,overlap);
y = y(1:length(sig(:,1)),:);