% Channel Esrimation
M = 10;
K = 3;
N = 16;
L = 4;
tau_p = 5;              % Pilot length
SNR_dB = 20;
SNR = 10.^(SNR_dB/10);
N0 = 10.^(-SNR_dB/10);
beta = 0;
teta = 0;
phi = 0;

% ===============================UE SIDE=====================
[ Xpf ] = genPilot( tau_p, K, N );

%  %======================Rayleigh Channel================
[Ht, Hf] = genChannel('Rayleigh',K,L,M,N,beta,teta,phi);

%=======================AWGN Noise======================
[ nt, nf ] = genAWGN ( M, tau_p, N );

Hft = zeros(M,K,N);
Hfy = zeros(M,K,N);
Hf_est = zeros(M,K,N);
for n = 1:N;  
    Hft(:,:,n) = Hf(:,:,n)';
    Hfy(:,:,n) = Hft(:,:,n)*Xpf(:,:,n);  
    Yp = Hfy + sqrt(N0)*nf;
    Hf_est(:,:,n) = Yp(:,:,n)/(Xpf(:,:,n));
end
error = abs(abs(estG)-abs(Hft));