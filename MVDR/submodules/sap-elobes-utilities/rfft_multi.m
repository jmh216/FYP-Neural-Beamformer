function x_f = rfft_multi(x,w_len,ov,ver)
% RFFT_MULTI gives the STFT of multichannel real data, based on the rfft function in [1].
%
% INPUTS:
% x:  	 	time-domain multichannel signal, with one column as one channel.
% w_len: 	window length 
% ov:		overlap
%
% OUTPUT:
% x_f:   	nMic*nFreq*nFrame STFT
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

if ver == 1508

	if ov == 0.75
	    w=hamming(w_len,'periodic')';     
	else
	    w=sqrt(hamming(w_len,'periodic'))'; 
	end

	% w=sqrt(hamming(w_len,'periodic'))'; 

	% w(end)=[]; % for now always use sqrt hamming window
	inc_len = fix(w_len*(1-ov));
	% w=w/sqrt(sum(w(1:inc_len:w_len).^2)); 
	w=0.5*hamming(w_len+1)/1.08; w(end)=[];
	[n_sig,nMic] = size(x);
	nFrame = max(fix((n_sig-w_len+inc_len)/inc_len),0);
	na=n_sig-w_len+inc_len-inc_len*nFrame+(nFrame==0)*(w_len-inc_len);       % number of samples left over
	nFrame = nFrame+(na>0);
	n_fft = w_len;
	y_f = zeros(nFrame,n_fft/2+1,nMic);

	for i_mic = 1:nMic
		y=enframe(x(:,i_mic),w,inc_len,'z');
		y_f(:,:,i_mic)=rfft(y,n_fft,2);
	end
else
	hlen = fix(w_len*(1-ov));
	[xlen,nMic] = size(x);
	rown = ceil((1+w_len)/2);            % calculate the total number of rows
	coln = 1+fix((xlen-w_len)/hlen);        % calculate the total number of columns
	y_f = zeros(coln,rown,nMic);           % form the stft matrix

	for i_mic = 1:nMic
		y_f(:,:,i_mic)=stft_z(x(:,i_mic), w_len, hlen, w_len).';
	end
end

if nMic>1
	x_f = permute(y_f,[3 2 1]);
else
	x_f = y_f.';
end