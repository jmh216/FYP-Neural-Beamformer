function [y,qout] = mkf_nn_smoothing_factor(sig,air,pp)
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
    [~,channel] = check_domain(pp,'channel',1);
    [~,fs] = check_domain(air,'fs',[]);
    
    if isreal(sig) %check whether the input is time-domain or stft-domain.
        sigSpec = rfft_multi(sig,wLen,overlap);
    else
        sigSpec = sig;
    end
    [~,nFreq,nFrame] = size(sigSpec);
    
    
    
    nOutChannel = length(channel);
    refChannel = channel(1);
    outSpec = zeros(nFreq,nFrame);
    bf_para = struct('beamformer','mvdr', 'win_len',wLen,'overlap',overlap);
    y_mvdr = beamforming(sig,air,bf_para);
    fs_nn = 10000;
    if fs_nn~=fs
        y = resample(y_mvdr(:,1),fs_nn,fs);
    end
    [ ~, binmasks, nfrm, nhopm ] = ComputeMaskEndtoEnd3({y}, 1);
    msk = covert_mask(y_mvdr(:,1),fs,fs_nn,binmasks{1},nfrm,nhopm,overlap,wLen);
    msk = msk';
    
    refSig = pp.ref_sig(:,1);
    lpcOrder = 2;

    sigRefSpec = rfft_multi(refSig,wLen,overlap);
    sigRefAmp = abs(sigRefSpec);
    ryyMat = stft2rxx(sigSpec,0.92);
    rvvMat = ryyMat;
    rvv = mean(ryyMat(:,:,:,1:50),4);

    len_m = 4;
    lpc_m = 8;
    refAmpTmp = sigRefAmp(:,1:lpc_m);
    [lpc_f,varw] = lpc(refAmpTmp',lpcOrder);
    A_lpc = zeros(lpcOrder,lpcOrder);
    A_lpc(2:end,1:end-1) = eye(lpcOrder-1);
    Pn = zeros(lpcOrder,lpcOrder,nFreq);
    PnTilde = zeros(lpcOrder,lpcOrder,nFreq);
    d = [1;zeros(lpcOrder-1,1)];
    xVec = permute(sigSpec(1,:,lpcOrder:-1:1),[3 2 1]);

    steerVec = rir2stv(air.h,nFreq);  

    for iFrame = 1:nFrame
        if mod(iFrame,len_m) == 0 && iFrame > lpc_m
            refAmpTmp = sigRefAmp(:,iFrame-lpc_m+1:iFrame);
            [lpc_f,varw] = lpc(refAmpTmp',lpcOrder);
        end

        sig_f = sigSpec(:,:,iFrame);
        ryy = ryyMat(:,:,:,iFrame);
        qq.ryy = ryy;
        qq.mask = msk(:,iFrame);
        [rvv,qq] = estnoicov_nn(sig_f,rvv,qq);
        rvvMat(:,:,:,iFrame) = rvv;

        for i_freq = 1:nFreq        
            A_lpc(1,:) = -lpc_f(i_freq,2:end);
            xVecAbs = A_lpc*abs(xVec(:,i_freq));%eq4
            if iFrame >= lpcOrder
                refPhase = sigRefSpec(i_freq,iFrame:-1:iFrame-lpcOrder+1).';
            else
                refPhase = [sigRefSpec(i_freq,iFrame:-1:1) zeros(1, lpcOrder-iFrame)].';
            end
            refPhase = refPhase./abs(refPhase+(1e-200));
            PhiK = diag(refPhase);
            xVecTmp = xVecAbs.*refPhase;%insert phase
            PnTildeTmp = A_lpc*PnTilde(:,:,i_freq)*A_lpc'+varw(i_freq)*d*d';%11
            PnTmp = PhiK*PnTildeTmp*PhiK';%12
            QMat = steerVec(:,i_freq)*[1 zeros(1,lpcOrder-1)];
            G = PnTmp*QMat'/(QMat*PnTmp*QMat'+rvv(:,:,i_freq));%9
            xVec(:,i_freq) = xVecTmp+G*(sig_f(:,i_freq)-QMat*xVecTmp);%6
            Pn(:,:,i_freq) = (eye(lpcOrder)-G*QMat)*PnTmp;%13
            PnTilde(:,:,i_freq) = PhiK'*Pn(:,:,i_freq)*PhiK;%14
            outSpec(i_freq,iFrame) = xVec(1,i_freq);
        end
    end

    outSpecMat = zeros(nFreq,nFrame,nOutChannel);
    
    for iChan = 1:nOutChannel
        RTF = steerVec(iChan,:)./steerVec(refChannel,:);
        outSpecMat(:,:,iChan) = outSpec.*repmat(RTF.',1,nFrame);
    end
    
    y = irfft_multi(permute(outSpecMat,[3,1,2]),wLen,overlap);
    y = y(1:length(sig(:,1)),:);
    
    qout.spp = msk;
    qout.rvv = rvvMat;
    qout.wLen = wLen;
    qout.overlap = overlap;

    
function bm = covert_mask(x,fs,fs_nn,bmm,nfrm,nhopm,overlap,wLen)
    of=1/(1-overlap); % overlap factor
    nf=wLen;      % fft length (samples)
    xtf=v_stft(x,nf,of); % apply STFT  
    
    
    bm=maskinterp(bmm,[fs/fs_nn*nfrm/nf 0 size(xtf,2)],[fs_nn/nhopm*nf/of/fs ((nf-1)*fs_nn/fs-nfrm+1)/nhopm/2+1 size(xtf,1)]);
    