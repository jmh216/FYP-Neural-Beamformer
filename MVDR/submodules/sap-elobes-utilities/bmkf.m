function out_sig = bmkf(mix_sig,para)
% BMKF uses the multichannel Kalman filter for speech enhancement in which the RTF is estimated blindly.
% 
% INPUTS:
% mix_sig: multichannel noisy input, one channel in one column.
% fs: sampling frequency
%
% OUTPUT:
% out_sig: single channel clean output signal
%
% 
% Wei Xue, Imperial College, May 03, 2017

if isstruct(para) || isempty(para)
    default_fields = {'win_len','overlap','lpc_w'};
    default_values = [128, 0.5, 1];
    n_fields = length(default_fields);
	for i_field = 1:n_fields
		if ~isfield(para,default_fields{i_field})
			eval([default_fields{i_field} '=' num2str(default_values(i_field)) ';']);
		else
			eval([default_fields{i_field} '=para.' default_fields{i_field} ';']);
		end
    end
end

sig_f_multi = rfft_multi(mix_sig,win_len,overlap);
[n_mic,n_freq,n_frame] = size(sig_f_multi);
ryy_mat = stft2rxx(sig_f_multi,0.92);
rvv = mean(ryy_mat(:,:,:,1:12),4);
spp = zeros(n_freq,1);
qq = struct('spp',spp,'phi_0',-2e3);
ref_sig = para.ref_sig(:,1);
lpc_order = 2;

sig_ref_f = rfft_multi(ref_sig,win_len,overlap);
sig_ref_mdl = abs(sig_ref_f);

len_m = 4;
lpc_m = 8;
ref_m_tmp = sig_ref_mdl(:,1:lpc_m);
[lpc_f,varw] = lpc(ref_m_tmp',lpc_order);
A_lpc = zeros(lpc_order,lpc_order);
A_lpc(2:end,1:end-1) = eye(lpc_order-1);
clean_f = zeros(2,n_freq,n_frame);
Pn = zeros(lpc_order,lpc_order,n_freq);
PnTilde = zeros(lpc_order,lpc_order,n_freq);
d = [1;zeros(lpc_order-1,1)];
xVec = permute(sig_f_multi(1,:,lpc_order:-1:1),[3 2 1]);

rtf = zeros(n_mic,n_freq);
ryy = ryy_mat(:,:,:,1);        
for i_freq = 1:n_freq
    rtf(:,i_freq) = power_eig(ryy(:,:,i_freq)-rvv(:,:,i_freq));
end

rtfqq = struct('alg',para.rtf_est,'rtfin',rtf,'qs',0.1);

for i_frame = 1:n_frame
    if mod(i_frame,len_m) == 0 && i_frame > lpc_m
        ref_m_tmp = sig_ref_mdl(:,i_frame-lpc_m+1:i_frame);
        [lpc_f,varw] = lpc(ref_m_tmp',lpc_order);
    end

    sig_f = sig_f_multi(:,:,i_frame);
    ryy = ryy_mat(:,:,:,i_frame);
    qq.ryy = ryy;
    [rvv,qq] = estnoicov_s(sig_f,rvv,qq);
    [rtfout,rtfqq] = estrtf(ryy,rvv,rtfqq);

	for i_freq = 1:n_freq        
        A_lpc(1,:) = -lpc_f(i_freq,2:end);
        A_lpc = A_lpc*lpc_w;
        xVecAbs = A_lpc*abs(xVec(:,i_freq));%eq4
        if i_frame >= lpc_order
            refPhase = sig_ref_f(i_freq,i_frame:-1:i_frame-lpc_order+1).';
        else
            refPhase = [sig_ref_f(i_freq,i_frame:-1:1) zeros(1, lpc_order-i_frame)].';
        end
        refPhase = refPhase./abs(refPhase+1e-200);
        PhiK = diag(refPhase);
        xVecTmp = xVecAbs.*refPhase;%insert phase
        PnTildeTmp = A_lpc*PnTilde(:,:,i_freq)*A_lpc'+varw(i_freq)*d*d';%11
        PnTmp = PhiK*PnTildeTmp*PhiK';%12
        QMat = rtfout(:,i_freq)*[1 zeros(1,lpc_order-1)];
        G = PnTmp*QMat'/(QMat*PnTmp*QMat'+rvv(:,:,i_freq));%9
        xVec(:,i_freq) = xVecTmp+G*(sig_f(:,i_freq)-QMat*xVecTmp);%6
        Pn(:,:,i_freq) = (eye(lpc_order)-G*QMat)*PnTmp;%13
        PnTilde(:,:,i_freq) = PhiK'*Pn(:,:,i_freq)*PhiK;%14
        clean_f(1,i_freq,i_frame) = xVec(1,i_freq);
        clean_f(2,i_freq,i_frame) = clean_f(1,i_freq,i_frame)*rtfout(2,i_freq)./rtfout(1,i_freq);
    end
end

clean_sig = irfft_multi(clean_f,win_len,overlap);
out_sig = clean_sig(1:length(mix_sig(:,1)),:);