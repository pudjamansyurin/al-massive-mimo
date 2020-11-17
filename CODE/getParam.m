function [ par ] = getParam( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% M = 50;               % Number of Tx antenna (in one BS)
par.M = 20:10:100;
par.N = 100;                % Number of occupied subcarrier
par.K = 10;                 % Number of Rx antenna (= number of UE)
% K = 10;
par.L = 4;                  % Channel tap frequency selective
par.beta = 1;
par.BPS = 2;                % (Bit/Symbol) Number of bits 
par.nBit = 2;               % Number bit per symbol
par.nCP = ceil(0.25*par.N);     % Number of cyclic Prefix (25% of NFFT)
% SNR_dB = 0:1:10;    % list of SNR [dB] values to be simulated
par.SNR_dB = 10;
par.Rn = 10.^(-par.SNR_dB/10);
par.SNR_L = 10^(par.SNR_dB/10);
par.FRM = 1;              % Number of data frame
par.tau_p = 20;
par.BPU = par.N*2;              % (Bit/User)  
par.NBPU = par.BPU*par.FRM;

par.QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];
par.symbol = par.QAM_symbol / sqrt(2); 
par.Code = { 
     'ZF' 
     'MRT' 
     'MMSE' 
};
par.Channel = { 
    'LOS' 
    'Rayleigh' 
};
end

