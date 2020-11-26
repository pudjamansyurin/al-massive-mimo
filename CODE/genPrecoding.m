function [A, C] = genPrecoding(type,K,M,Hf,S,N,Pc,SNR_L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


factor = 0;
if strcmp(type, 'MMSE')
    % MMSE Precoding
    A =  Hf'/(Hf*Hf'+ (1/Pc)*diag(ones(K,1)));  % Precoding matrix 
elseif strcmp(type, 'ZF')   
    A = Hf'/(Hf*Hf');
else
    % MRT Precodng
    A = 1/M*Hf';
end
% factor = factor+trace(A*A')/N;    %Precoding Factor
% A = A./sqrt(factor); % Scaled output
% Precoded Vector
C = A*S;  
end

