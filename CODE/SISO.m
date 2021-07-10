%%%SISO
M = 1;              % Number of occupied subcarrier
K = 1;
BPS = 2;                % (Bit/Symbol) Number of bits 
nBit = 2;               % Numer \bit per symbol
SNR_dB = 10;    % list of SNR [dB] values to be simulated
SNR_L = 10^(SNR_dB(length(SNR_dB))/10);
BPU = N*2;              % (Bit/User)  
jBit = 100;
QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];
symbol = QAM_symbol / sqrt(2);
B = randi([0, 1], [K jBit]);
Tx = modQAM(BPS, symbol, B);
Rx = awgn(Tx, SNR_L);
slot = reshape(Rx, 2, []);
demodulasi = demodQAM(slot, symbol)