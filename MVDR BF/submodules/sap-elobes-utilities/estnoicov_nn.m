function [rvv,qq] = estnoicov_nn(sig_f,rvv,qq,pp)
    % ESTNOICOV_S returns the souden's method [1] to estimate the
    % frequency-domain noise covariance matrix.
    %
    % INPUTS:
    % sig_f: nMic*n_fft frequency-domain multichannel signal vector of one frame
    % rvv:	 nMic*nMic*n_fft noise covariance matrix
    % qq:	 variables that need to be updated frame-wise
    % 		 qq.ryy: signal covariance matrix
    % 		 qq.spp: speech presence probability of last frame
    % 
    % pp:	 parameters used in [1].
    %		 pp.epsi: significance level in (12)
    % 		 pp.L:	 the number of periodograms for average in (13)
    %		 pp.K1:	 defines the size of normalized Hann window in (19)
    % 		 pp.alpha_p: smoothing parameter for multichannel speech presence probability
    % 		 pp.alpha_v: smoothing parameter in (25)
    %		
    % OUTPUT:
    % Rvv:  updated nMic*nMic*n_fft noise covariance matrix
    % spp:	speech presence probability of one frame
    %
    % REFERENCES:
    % [1] Souden, M., Jingdong Chen, J. Benesty, and S. Affes. an Integrated Solution for Online Multichannel Noise Tracking and Reduction.
    % IEEE Transactions on Audio, Speech, and Language Processing 19, no. 7 (September 2011).
    % 
    % Wei Xue, Imperial College, Feb 14, 2017
    if nargin<=3
        pp = [];
    end
    [~,alpha_v] = check_domain(pp,'alpha_v',0.92);
    msk = qq.mask;
    [~, n_fft] = size(sig_f);
    alpha_v_hat = alpha_v + (1- alpha_v)*msk;
    
    
    for i_freq = 1:n_fft
        rvv(:,:,i_freq) = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f(:,i_freq)*sig_f(:,i_freq)'); %eq24
    end