function [ SINR,Sig,I,Noise,Pt,AA,Sig2] = calcSINR (A, Hf, SNR_dB, K, Pc, SNR_L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
AA = A;
Pt = zeros(1,K);
for Ki = 1:K;
   Pt(Ki) = Pc*mean(abs(A(:,Ki)).^2);
end

% Signal Power
Sig =   SNR_L*diag(abs(Hf*A).^2);
Sig2 = abs(Hf*A).^2;
% Interference Power
I = SNR_L*sum(abs(Hf*A).^2,2)-Sig; 
% Noise Power
Noise = repmat(10.^(-SNR_dB/10),K,1);
% SINR of ZF at all user (simulation)
SINR = Sig./(I+Noise);