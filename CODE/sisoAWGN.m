clear all
clc;
SNR = 0:1:15;
SNR_L = 10.^(SNR/10);
M = 4;
x = 3*log2(M)/(M-1)*SNR_L;
mimo = [0.25 0.23 0.2 0.15 0.13 0.1 0.08 0.06 0.04 0.025 0.015 0.006 0.003 0.001 0.0003 0.00003]
sisoawgn = 4*qfunc(sqrt(x)); % BER SISO AWGN
% semilogy(SNR,sisoawgn,'--');
% % hold on;
semilogy(SNR,mimo);
legend('SISO AWGN')
