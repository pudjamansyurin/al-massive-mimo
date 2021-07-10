function [ SINR,Sig,I,Noise] = calcSINR (A, Hf, SNR_dB, K, Pc, SNR_L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Signal Power
Sig = diag(abs(Hf*A).^2);
% Interference Power
I = sum(abs(Hf*A).^2,2)-Sig; 
% Noise Power
Noise = repmat(10.^(-SNR_dB/10),K,1);
% SINR of ZF at all user (simulation)
SINR = Sig./(I+Noise);