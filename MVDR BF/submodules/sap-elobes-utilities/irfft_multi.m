function x = irfft_multi(x_f,w_len,ov,ver)
% IRFFT_MULTI gives the inverse DFT of multichannel data.
%
% INPUTS:
% x_f:  	nMic*n_fft*nFrame DFT
% w_len: 	window length 
% ov:		overlap
%
% OUTPUT:
% x:   	 	n_sig*nMic time-domain multichannel signal
%
% REFERENCES: 
% [1] D. M. Brookes, VOICEBOX: A speech processing toolbox for MATLAB, http://www.ee.ic.ac.uk/hp/staff/dmb/
% voicebox/voicebox.html, 1997-2017.

% Wei Xue, Imperial College, Feb 14, 2017

if nargin<4
	ver = 1508;
else
	ver = 1808;
end

if ndims(x_f)<3
	x_f = permute(x_f,[3 1 2]);
end
[nMic,n_fft,nFrame] = size(x_f);

if ver == 1508

	if ov == 0.75
	    w=hamming(w_len,'periodic')';     
	else
	    w=sqrt(hamming(w_len,'periodic'))'; 
	end
	w=sqrt(hamming(w_len,'periodic'))'; 

	% w(end)=[]; % for now always use sqrt hamming window
	inc_len = fix(w_len*(1-ov));
	w = w/sqrt(sum(w(1:inc_len:w_len).^2));

	
	n_sig = nFrame*inc_len-inc_len+w_len;
	x = zeros(n_sig,nMic);

	for i_mic = 1:nMic
		x_f_mic = permute(x_f(i_mic,:,:),[3 2 1]);
		x_mic = overlapadd(irfft(x_f_mic,w_len,2),w,inc_len);
		x(:,i_mic) = x_mic(:);
	end
else
	coln = size(x_f, 3);
	hlen = fix(w_len*(1-ov));
	nfft = w_len;
	n_sig = nfft + (coln-1)*hlen;
	x = zeros(n_sig,nMic);
	for i_mic = 1:nMic
		x_f_mic = permute(x_f(i_mic,:,:),[3 2 1]);
		x_mic = istft_z(x_f_mic.', hlen, nfft);
		x(:,i_mic) = x_mic(:);
	end
end