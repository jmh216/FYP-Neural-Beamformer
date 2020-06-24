function [y,qout] = beamforming(sig,air,pp)
% BEAMFORMING denoises signal using frequency-domain beamformers.
%
% INPUTS:
% sig:  time-domain multichannel signal, with one column as one channel.
% air:	acoutic environment
%       air.fs: sampling rate
%       air.h: room impulse responses(RIRs). The RIRs are only used in[3]
%       to compute the relative transfer function of early RIRs.
%       OPTIONAL air.noi: if air contains this field, the algorithm use
%       oracle noise information for mvdr/mwf.
%       Other fields of air may change depending on whether the scenario is
%       simulated or generated using the Oldenburg database.
% pp:   parameters of the algorithm. Default values are shown in brackets.
%       pp.wLen(256): Window size of short-time-Fourier-transform
%       (STFT). 
%       pp.overlap(0.5): Overlap between frames in range [0,1]
%       pp.f_debug(0):   Whether operates in the debug mode 
%       pp.f_rtf(0):     The steering vectors are computed using RIRs if f_rtf == 1.
%                        Otherwise, use the geometry information contained in "air" instead. 
%       pp.beamformer('mwf'): choose the beamformer type
%       pp.noi_est('s'): If 's', use souden's method [4] for noise
%       covariance matrix estimation. If 'h', use [3] instead.
%       pp.channel(1): the reference channel
% OUTPUTS:
% y:    time-domain single-channel beamformer output
% w_mat:nMic*nFreq*nFrame matrix consisting of frequency-domain filter
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
if nargin<3
    pp = [];
end

[~,wLen] = check_domain(pp,'wLen',256);
[~,overlap] = check_domain(pp,'overlap',0.5);
[~,smoothFactor] = check_domain(pp,'smoothFactor',0.8);
[~,channel] = check_domain(pp,'channel',1);
[~,bfName] = check_domain(pp,'beamformer','mwf');

if isreal(sig) %check whether the input is time-domain or stft-domain.
    sigSpec = rfft_multi(sig,wLen,overlap);
else
    sigSpec = sig;
end
[~,nFreq,nFrame] = size(sigSpec);

if is_instring(bfName,'mwf')    
    ppbf = struct('channel',channel,'w',1);
elseif is_instring(bfName,'mvdr') 
    steerVec = rir2stv(air.h,nFreq);           
end

%Initialize ryy and rvv, set parameters
ryyMat = stft2rxx(sigSpec,smoothFactor);

if isfield(air,'noi') %non-blind, the noise signal is a field of pp
    noi = air.noi;
    if isreal(noi) %check whether the input is time-domain or stft-domain.
        noiSpec = rfft_multi(noi,wLen,overlap);
    else
        noiSpec = noi;
    end
	rvvMat = stft2rxx(noiSpec,smoothFactor);
	noi_flg = 1;
else %blind
    rvv = mean(ryyMat(:,:,:,1:50),4);
    qq = struct('spp',zeros(nFreq,1),'phi_0',-2e3);
    rvvMat = ryyMat;
	noi_flg = 0;
end

nOutChannel = length(channel);
refChannel = channel(1);
outSpec = zeros(nFreq,nFrame);
wfGain = zeros(nFreq,nFrame);
minVal = 10^(-25/10);

sppMat = zeros(nFreq,nFrame);
%framewise processing
for iFrame = 1:nFrame
	sigFFT = sigSpec(:,:,iFrame);
    ryy = ryyMat(:,:,:,iFrame);
	%_______________Noise estimation_______________
    if noi_flg
        rvv = rvvMat(:,:,:,iFrame); %use oracle noise       
    else
        qq.ryy = ryy;
        [rvv,qq] = estnoicov(sigFFT,rvv,qq);     
        sppMat(:,iFrame) = qq.spp;
        rvvMat(:,:,:,iFrame) = rvv;
    end

    %_______________beamforming_______________
    for iFreq = 1:nFreq
        switch bfName %computing beamformer coefficients, default: mwf
            case {'mvdr','mvdr-wf'}
                h_bf = h_mvdr(rvv(:,:,iFreq),steerVec(:,iFreq)./steerVec(refChannel,iFreq)); %mvdr
            case 'sp-mvdr'
                h_bf = h_sp_mvdr(ryy(:,:,iFreq),rvv(:,:,iFreq),steerVec(:,iFreq)./steerVec(refChannel,iFreq)); %mvdr
            case 'sp-mwf'
                h_bf = h_sp_mwf(ryy(:,:,iFreq),rvv(:,:,iFreq),ppbf); %wmf
            otherwise
                h_bf = h_mwf(ryy(:,:,iFreq),rvv(:,:,iFreq),ppbf); %wmf
        end

        outSpec(iFreq,iFrame) = h_bf'*sigFFT(:,iFreq); %spatial filtering

        if is_instring(bfName,'mvdr-wf')
            pyyo = h_bf'* ryy(:,:,iFreq)*h_bf;
            pvvo = h_bf'* rvv(:,:,iFreq)*h_bf;
            wfGainTF = max(minVal,1-real(pvvo/pyyo));
            outSpec(iFreq,iFrame) = outSpec(iFreq,iFrame).*wfGainTF;
            wfGain(iFreq,iFrame) = wfGainTF;
        end
    end    
end

outSpecMat = zeros(nFreq,nFrame,nOutChannel);

for iChan = 1:nOutChannel
    RTF = steerVec(iChan,:)./steerVec(refChannel,:);
    outSpecMat(:,:,iChan) = outSpec.*repmat(RTF.',1,nFrame);
end

y = irfft_multi(permute(outSpecMat,[3,1,2]),wLen,overlap);
y = y(1:length(sig(:,1)),:);

qout.spp = sppMat;
qout.channel = channel;
qout.rvv = rvvMat;
qout.wLen = wLen;
qout.overlap = overlap;
qout.wfGain = wfGain;