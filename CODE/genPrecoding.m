function [A, C, factor] = genPrecoding(type,K,M,Hf,S,N,Pc,SNRo)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

SNRd = 10^(SNRo/10);
factor = 0;
if strcmp(type, 'MMSE')
    % MMSE Precoding
    A =  Hf'/(Hf*Hf'+ (1/SNRd)*diag(ones(K,1)));  % Precoding matrix 
elseif strcmp(type, 'ZF')   
    A = Hf'/(Hf*Hf');
else
    % MRT Precodng
    A = 1/M*Hf';
end

factor = mean(trace(A*A'));    %Precoding Factor
% factor = factor + trace(A*A')/N;
A = A./sqrt(factor); % Scaled output
% % Precoded Vector
C = A*S;  
end