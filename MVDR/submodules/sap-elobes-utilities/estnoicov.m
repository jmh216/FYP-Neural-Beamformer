function [rvv,qq] = estnoicov(sig_f,rvv,qq,pp)

%--- Check and unpack parameters ---
if nargin<4
    pp = [];
end

[~,cmode] = check_domain(pp,'cmode','s');

switch cmode
    case 's'
        [rvv,qq] = estnoicov_s(sig_f,rvv,qq,pp);
    case 'h'
        [rvv,qq] = estnoicov_h(sig_f,rvv,qq);
    case 't'
        [rvv,qq] = estnoicov_t(sig_f,rvv,qq,pp);        
end


function [rvv,qq] = estnoicov_t(sig_f,rvv,qq,pp)
% ESTNOICOV_S returns the Maja Taseska's method [1] to estimate the
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
% [1] Taseska, Maja, and Emanuel A. P. Habets. â€œMMSE-Based Blind Source Extraction
% in Diffuse Noise Fields Using a Complex Coherence-Based a Priori SAP Estimator.ï¿?? 
% In Acoustic Signal Enhancement; Proceedings of IWAENC 2012; International Workshop on, 1â€?4, 2012.

% 
% Wei Xue, Imperial College, Mar 10, 2017

%--- Check and unpack parameters ---
if nargin<4
    pp = [];
end

if isstruct(pp) || isempty(pp)
    default_fields = {'alpha_p','alpha_v'};
    default_values = [0.6, 0.92];%default values
    n_fields = length(default_fields);
	for i_field = 1:n_fields
		if ~isfield(pp,default_fields{i_field})
			eval([default_fields{i_field} '=' num2str(default_values(i_field)) ';']);
		else
			eval([default_fields{i_field} '=pp.' default_fields{i_field} ';']);
		end
	end
else
	error('please specify the parameters.');	
end

ryy = qq.ryy;
spp = qq.spp;
noi_psd_mat = qq.noi_mat;
X_psda = qq.sig_psda;
X_psdb = qq.sig_psdb;
X_cpsd = qq.sig_cpsd;
d = qq.d;
fs = qq.fs;
c = 343;

%--- Other default parameters ---
[nMic,nFreq] = size(sig_f);
n_fft = (nFreq-1)*2;
max_val = 1e4;
min_sap = 0.4;
l_min = 0.2;
l_max = 0.8;
rho = 2;
c_val = 3;
tmp_val = 10^(c_val*rho/10);
tmp_val_w = 10^(0*rho/10);
wave_number = 2*pi*((0:nFreq-1)'/n_fft*fs)/c;
r_diff = sin(wave_number*d)./(wave_number*d+eps);

%--- Allocate memory ----
phi_hat = zeros(nFreq,1);
eta_hat = zeros(nFreq,1);
beta_hat = zeros(nFreq,1);

%--- Start Processing ----
for i_freq = 1:nFreq
	sig_f_vec = sig_f(:,i_freq);	
	ryy_i_freq = ryy(:,:,i_freq);
	rvv_i_freq = rvv(:,:,i_freq);

	rxx_i_freq = ryy_i_freq-rvv_i_freq;%2.2.a
	inv_rvv = inv(rvv_i_freq);
	phi_hat(i_freq) = real(trace(inv_rvv*ryy_i_freq));%2.2.c
	eta_hat(i_freq) = max(phi_hat(i_freq)-nMic,0);% 2.2.d
	beta_hat(i_freq) = real(sig_f_vec'*inv_rvv*rxx_i_freq*inv_rvv*sig_f_vec);%2.2.e
end		


S_psda = max(X_psda-noi_psd_mat(:,1),eps);
S_psdb = max(X_psdb-noi_psd_mat(:,2),eps);
complex_coherence = X_cpsd./((S_psda+S_psdb)/2);
mu_dir = angle(X_cpsd);%eq 14
srr = max(real((r_diff-complex_coherence)./(complex_coherence - exp(1i*mu_dir))),eps);
sap = l_min+(l_max-l_min).*tmp_val./(tmp_val+srr.^rho);


%__________________________________MC-SPP and noise estimation in the first and second stage(2.4-3.2)________________________________________
spp_0 = 1./(1+(sap./(1-sap)).*(1+ eta_hat).* min(real(exp(-beta_hat./(1+eta_hat))),max_val));%using val_max to avoid inf
spp_1 = alpha_p*spp + (1-alpha_p)*spp_0;
alpha_v_hat = alpha_v + (1- alpha_v)*spp_1;

for i_freq = 1:nFreq
	rvv(:,:,i_freq) = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f(:,i_freq)*sig_f(:,i_freq)'); %eq24
end

%--- Re-package variables---
qq.spp = spp;
qq.w = 5*tmp_val_w./(tmp_val_w+srr.^rho);


function [rvv,qq] = estnoicov_h(sig_f,rvv,qq)
% ESTNOICOV_H returns the Hendriks's method to estimate the frequency-domain noise covariance matrix.
%
% INPUTS:
% sig_f: nMic*n_fft frequency-domain multichannel signal vector of one frame
% rvv:	 nMic*nMic*n_fft noise covariance matrix
% qq:	 variables that need to be updated frame-wise
% 		 qq.ryy: signal covariance matrix
% 
% pp:	 parameters used in [1].
%		
% OUTPUT:
% Rvv:   updated nMic*nMic*n_fft noise covariance matrix
% spp:	 speech presence probability of one frame
%
% REFERENCES:
% [1] Hendriks, R.C., and T. Gerkmann. Noise Correlation Matrix Estimation for
% Multi-Microphone Speech Enhancement. IEEE Transactions on Audio, Speech, and Language 
% Processing 20, no. 1 (2012): 223-233. doi:10.1109/TASL.2011.2159711.

% 
% Wei Xue, Imperial College, Feb 16, 2017

%--- Check and unpack parameters ---
if ~isfield(qq,'eq')
    qq.eq = 11;
end

rpy = qq.rpy;
noi_psd_mat = qq.noi_mat;
steerVec = qq.steerVec;

[nMic,n_fft] = size(sig_f);

%--- Estimate noise psd of each channel ---
for n = 1:nMic
	rvv(n,n,:) = noi_psd_mat(n,:);
end

%--- Start Processing ----
for i_freq = 1:n_fft
	sig_f_vec = sig_f(:,i_freq);
	steering_vec = steerVec(:,i_freq);

    if 10*log10(abs(sig_f_vec(1))^2/rvv(1,1,i_freq))>7.8
        alpha_p = 0.01;
    else
        alpha_p = 0.1;
    end
			
	for n = 1:nMic
		for m = n+1:nMic
			dnm = steering_vec(n)/steering_vec(m);
			pnm = sig_f_vec(n)-dnm*sig_f_vec(m);
			rpy(n,m,i_freq) = (1-alpha_p)*rpy(n,m,i_freq)+alpha_p*pnm*conj(sig_f_vec(m));
			
			dmn = steering_vec(m)/steering_vec(n);
			pmn = sig_f_vec(m)-dmn*sig_f_vec(n);
			rpy(m,n,i_freq) = (1-alpha_p)*rpy(m,n,i_freq)+alpha_p*pmn*conj(sig_f_vec(n));

			if qq.eq == 11
				enm = rpy(n,m,i_freq)+dnm*rvv(m,m,i_freq);
				emn = rpy(m,n,i_freq)+dmn*rvv(n,n,i_freq);
				rvv(n,m,i_freq) = (enm+conj(emn))/2;
				rvv(m,n,i_freq) = conj(rvv(n,m,i_freq));
			else
				rvv(n,m,i_freq) = (rpy(n,m,i_freq)+rpy(m,n,i_freq))/2;
				rvv(m,n,i_freq) = conj(rvv(n,m,i_freq));
			end
		end
    end
end

if qq.eq ~= 11
    for n = 1:nMic
        rvv(n,n,:) = 0;
    end
end

%--- Re-package variables---
qq.rpy = rpy;
    
function [rvv,qq] = estnoicov_s(sig_f,rvv,qq,pp)
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

[~,epsi] = check_domain(pp,'epsi',0.01);
[~,L] = check_domain(pp,'L',32);
[~,K1] = check_domain(pp,'K1',15);
[~,alpha_p] = check_domain(pp,'alpha_p',0.6);
[~,alpha_v] = check_domain(pp,'alpha_v',0.92);

ryy = qq.ryy;
spp = qq.spp;

%--- Other default parameters ---
max_val = 1e4;
min_sap = 0.4;

%--- Allocate memory ----
[nMic, n_fft] = size(sig_f);
phi_hat = zeros(n_fft,1);
eta_hat = zeros(n_fft,1);
beta_hat = zeros(n_fft,1);
phi_inst = zeros(n_fft,1);
sap_local = zeros(n_fft,1);

%---compute phi_0 and phi_0_h in the first frame---
if qq.phi_0 < 0
    [phi_0, phi_0_h] = solve_cum(nMic,L,epsi);
    qq.phi_0 = phi_0;
    qq.phi_0_h = phi_0_h;
else
    phi_0 = qq.phi_0;
    phi_0_h = qq.phi_0_h;
end

w_global = hann(K1);
w_global = w_global./sum(w_global);

%--- Start Processing ----
for i_freq = 1:n_fft
	sig_f_vec = sig_f(:,i_freq);	
	ryy_i_freq = ryy(:,:,i_freq);
	rvv_i_freq = rvv(:,:,i_freq);

	rxx_i_freq = ryy_i_freq-rvv_i_freq;%2.2.a
	inv_rvv = inv(rvv_i_freq);
	phi_inst(i_freq) = real(sig_f_vec'*inv_rvv*sig_f_vec);%2.2.b
	phi_hat(i_freq) = real(trace(inv_rvv*ryy_i_freq));%2.2.c
	eta_hat(i_freq) = max(phi_hat(i_freq)-nMic,0);% 2.2.d
	beta_hat(i_freq) = real(sig_f_vec'*inv_rvv*rxx_i_freq*inv_rvv*sig_f_vec);%2.2.e

	%---eq18---
	if phi_inst(i_freq) < phi_0 && phi_hat(i_freq) < nMic
		sap_local(i_freq) = 1;
	elseif phi_inst(i_freq) < phi_0 && phi_hat(i_freq) >= nMic && phi_hat(i_freq)<phi_0_h
		sap_local(i_freq) = (phi_0_h-phi_hat(i_freq))/(phi_0_h-nMic);
	else
		sap_local(i_freq) = 0;
	end
end		

phi_global = conv(phi_inst,w_global,'same'); %eq19
phi_frame = mean(phi_inst); %eq20
sap_global = phi_global < phi_0;
sap_frame = phi_frame < phi_0;
sap = min(sap_local.*sap_global.*sap_frame,0.99);
sap = max(sap,min_sap);

%__________________________________MC-SPP and noise estimation in the first and second stage(2.4-3.2)________________________________________
spp_0 = 1./(1+(sap./(1-sap)).*(1+ eta_hat).* min(real(exp(-beta_hat./(1+eta_hat))),max_val));%using val_max to avoid inf
spp_1 = alpha_p*spp + (1-alpha_p)*spp_0;
alpha_v_hat = alpha_v + (1- alpha_v)*spp_1;


 %________________Second Stage_____________________________
for i_freq = 1:n_fft
	sig_f_vec = sig_f(:,i_freq);
	rvv_i_freq = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f_vec*sig_f_vec');
	ryy_i_freq = ryy(:,:,i_freq);
	
	rxx_i_freq = ryy_i_freq - rvv_i_freq; %2.2.a
	inv_rvv = inv(rvv_i_freq);
	phi_inst(i_freq) = real(sig_f_vec'*inv_rvv*sig_f_vec);%2.2.b
	phi_hat(i_freq) = real(trace(inv_rvv*ryy_i_freq));%2.2.c
	eta_hat(i_freq) = max(phi_hat(i_freq)-nMic,0);% 2.2.d
	beta_hat(i_freq) = real(sig_f_vec'*inv_rvv*rxx_i_freq*inv_rvv*sig_f_vec);%2.2.e

	%---eq18---
	if phi_inst(i_freq) < phi_0 && phi_hat(i_freq) < nMic
		sap_local(i_freq) = 1;
	elseif phi_inst(i_freq) < phi_0 && phi_hat(i_freq) >= nMic && phi_hat(i_freq)<phi_0_h
		sap_local(i_freq) = (phi_0_h-phi_hat(i_freq))/(phi_0_h-nMic);
	else
		sap_local(i_freq) = 0;
	end
end

phi_global = conv(phi_inst,w_global,'same'); %eq19
phi_frame = mean(phi_inst); %eq20
sap_global = phi_global < phi_0;
sap_frame = phi_frame < phi_0;
sap = min(sap_local.*sap_global.*sap_frame,0.99);
sap = max(sap,min_sap);

%---Final Estimation---
spp = 1./(1+(sap./(1-sap)).*(1+ eta_hat).* min(real(exp(-beta_hat./(1+eta_hat))),max_val));
alpha_v_hat = alpha_v + (1- alpha_v)*spp;

for i_freq = 1:n_fft
	rvv(:,:,i_freq) = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f(:,i_freq)*sig_f(:,i_freq)'); %eq24
end

%--- Re-package variables---
qq.spp = spp;

function [y0, y1] = solve_cum(nMic,L,eps0)

F_phi = @(x) (x/L).^nMic*L*gamma(L)/(gamma(nMic+1)*gamma(L-nMic+1)).*hypergeom([nMic,L+1],nMic+1,-x/L).*double(x>=1); %eq15
a_F = 2*nMic*L;
B_F = (L+L-2*nMic-1)*(L-1)/((L-2*nMic-3)*(L-2*nMic));
b_F = 4+(a_F+2)/(B_F-1);
c_F = a_F*(b_F-2)/(b_F*(L-2*nMic-1));
F_phi_h = @(x) 1.*betainc(a_F*x./(a_F*x+b_F),a_F/2,b_F/2).*double(x>=1);

y0 = solve_f(F_phi,1-eps0,1e-5);
y1 = c_F*solve_f(F_phi_h,1-eps0,1e-5);