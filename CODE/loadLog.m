load 'logs/20200930T123035.mat'

semilogy(SNR_dB, ber.ratio,'-rs')
hold on;
semilogy(SNR_dB, ber.ratio_MRT,'-bs')
title('Bit Error Rate MU-Massive MIMO');
xlabel('SNR(dB)');
ylabel ('BER');
legend('ZF','MRT');
grid on;