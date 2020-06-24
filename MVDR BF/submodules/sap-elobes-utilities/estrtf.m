function [rtfout,qq] = estrtf(ryy,rvv,qq)
% ESTRTF returns the estimated relative transfer function (RTF) using the covariance whitening method.
%
% INPUTS:
% ryy: 		n_mic*n_mic*n_fft frequency-domain multichannel signal vector of one frame
% rvv:	 	n_mic*n_mic*n_fft noise covariance matrix
% qq:		optional parameters
% rtfin:	n_mic*n_fft rtf matrix of the previous frame
% um:		n_mic*n_fft matrix containing the generalized eigenvectors
% OUTPUT:
% rtfout:	n_mic*n_fft rtf matrix of the current frame
%
% REFERENCES:
% [1] Varzandeh, Reza, Maja Taseska, and E. Habets. "An iterative multichannel subspace-based covariance
% subtraction method for relative transfer function estimation." Hands-free Speech Communications and Microphone Arrays (HSCMA), 2017. IEEE, 2017.
% 
% Wei Xue, Imperial College, April 25, 2018

alg = qq.alg;
[n_mic,~,n_freq] = size(ryy);

if strcmp(lower(alg),'cw')
	um = qq.um;
	for i_freq = 1:n_freq
		ryyi = ryy(:,:,i_freq);
		rvvi = rvv(:,:,i_freq);
		umax = rvvi\ryyi*um(:,i_freq); %eq.10
		umax = umax/norm(umax); %eq.10
	    um(:,i_freq) = umax;
		tmp = rvvi*umax; %eq.9
		rtfout(:,i_freq) = tmp/tmp(1); %eq.9
	end
	qq.um = um;

elseif strcmp(lower(alg),'cs-evd')
	rtfin = qq.rtfin;
	for i_freq = 1:n_freq
		ryyi = ryy(:,:,i_freq);
		rvvi = rvv(:,:,i_freq);
		rxxi = ryyi - rvvi;
		[u,v] = eig(rxxi);
		[~,indx] = max(abs(diag(v)));
		if min(real(diag(v)))<0
			rtfout(:,i_freq) = rtfin(:,i_freq);
		else
			rtfout(:,i_freq) = u(:,indx)/u(1,indx); %sec 3.3
		end
	end
	qq.rtfin = rtfout;
elseif strcmp(lower(alg),'cs')
	rtfin = qq.rtfin;
	for i_freq = 1:n_freq
		ryyi = ryy(:,:,i_freq);
		rvvi = rvv(:,:,i_freq);
		rxxi = ryyi - rvvi;
		if min(real(diag(rxx)))<0
			rtfout(:,i_freq) = rtfin(:,i_freq);
		else
			rtfout(:,i_freq) = rxxi(:,1)/rxxi(1,1); %eq.6
		end
	end	
	qq.rtfin = rtfout;
end