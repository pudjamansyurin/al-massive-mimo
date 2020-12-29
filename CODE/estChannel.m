function [ Hf_est ] = estChannel( Hf, Xpf, N0, nf, tau_p, SNR_L )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

Hft = Hf';
Hfy = Hft*Xpf;
Yp = Hfy + sqrt(N0)*nf;
Yp_aksen = Yp*Xpf;
% Hf_est  = sqrt(tau_p*SNR_L)/(1+tau_p*SNR_L)*Yp_aksen;
Hf_est = Yp/Xpf;
Hf_est = Hf_est';

end

