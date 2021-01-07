function [ Hf_est ] = estChannel( Hf, Xpf, N0, nf, tau_p, SNRo,fi2)

SNRul = SNRo;
SNRuldB = 10^(SNRul/10) ;
Hft = Hf';
Hfy = sqrt(SNRuldB)*Hft*Xpf;
Yp = Hfy + sqrt(N0)*nf;
Yp_aksen = Yp*fi2;
Hf_est  = (sqrt(tau_p*SNRuldB)/(1+tau_p*SNRuldB))*Yp_aksen;
% Hf_est = Yp/Xpf;
Hf_est = Hf_est';

end

