function Rxx = stft2rxx(X,alpha_s,mode)
% STFT2RXX returns the spatial covariance matrix using exponential
% smoothing.
%
% INPUTS:
% X: 3-D nMic*n_freq*nFrame STFT matrix.
% alpha_s: exponential smoothing factor.
%
% OUTPUT:
% Rxx:  4-D nMic*nMic*n_freq*nFrame spatial covariance matrix
%
% Wei Xue, Imperial College, Feb 14, 2017

if nargin<3
    mode = 's';
end

[nMic,n_freq,nFrame] = size(X);


if strcmp(mode,'s')
    Rxx = zeros(nMic,nMic,n_freq,nFrame);

    n_init = 12;
    Rxx_init = zeros(nMic,nMic,n_freq);
    for i_frame = 1:n_init
        X_frame = X(:,:,i_frame);
        for i_freq = 1:n_freq
            Rxx_init(:,:,i_freq) = Rxx_init(:,:,i_freq)+X_frame(:,i_freq)*X_frame(:,i_freq)'; 
        end
    end

    for i_freq = 1:n_freq
        Rxx(:,:,i_freq,1) = Rxx_init(:,:,i_freq)/n_init;
    end

    %--- exponential smoothing ---
    for i_frame = 2:nFrame
        X_frame = X(:,:,i_frame);
        for i_freq = 1:n_freq
            Rxx(:,:,i_freq,i_frame) = alpha_s*Rxx(:,:,i_freq,i_frame-1) + (1-alpha_s)*(X_frame(:,i_freq)*X_frame(:,i_freq)');
        end
    end
else
    n_init = nFrame;
    Rxx_init = zeros(nMic,nMic,n_freq);
    Rxx = zeros(nMic,nMic,n_freq);
    for i_frame = 1:n_init
        X_frame = X(:,:,i_frame);
        for i_freq = 1:n_freq
            Rxx_init(:,:,i_freq) = Rxx_init(:,:,i_freq)+X_frame(:,i_freq)*X_frame(:,i_freq)'; 
        end
    end

    for i_freq = 1:n_freq
        Rxx(:,:,i_freq) = Rxx_init(:,:,i_freq)/n_init;
    end
end