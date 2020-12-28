function [ Hf_est ] = estChannel( Hf, Xpf, N0, nf )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

Hft = Hf';
Hfy = Hft*Xpf;
Yp = Hfy + sqrt(N0)*nf;
Yp = Yp*conj(Xpf);
Hf_est = Yp/Xpf;
% Hf_est = Yp/20;
Hf_est = Hf_est';

end

