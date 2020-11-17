% Generate Haar distributed random unitary matrices
% From book Acta Numerica

clear all; clc;         % Clear screen
tau_p = 5;              % Pilot length
M = 10;                 % Number of Tx antenna (in one BS)
K = 3;                  % Number of Rx antenna (= number of UE)
L = 4;                  % Channel tap frequency selective
v = 1 / sqrt(L);        % Variance channel for frequency selective
vn =1;                 % Variance noise
mu = 0;                 % Mean of noise
beta = 1;               % Large scale fadding Coefficient Of Rayleigh Channel
SNR_ul = 128;           % Uplink SNR (dB)
SNR_dl = SNR_ul;        % Download SNR (dB)
nBit = 2;               % Number of bits per chunk
QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];

% ===================== UE SIDE ==========================
% GENERATE PILOT -----------------------------------------
% Generate a random complex matrix 
X = complex(rand(tau_p), rand(tau_p)) / sqrt(2);
% Factorize the matrix
[Q, R] = qr(X);
% Unitary matrix M x N
fi = Q(:, 1:K);
% Verification of unitary matrix
verify = ctranspose(fi) * fi;
% Generate the pilot
Xp = sqrt(tau_p) * ctranspose(fi);


% ===================== UE > BS ==========================
% Standard deviation
sigma = sqrt(v);        
% Generate the Rayleigh channel from UE to BS
H = sigma .* complex(rand(M, K), rand(M, K)) + mu;
% Add large scale fadding to the Rayleigh channel
G = sqrt(beta) * H;
% Generate noise (AWGN)
Wp = sqrt(vn) .* complex(rand(M, tau_p), rand(M, tau_p)) + mu;


% ===================== BS SIDE ==========================
% Received pilot signal
Yp = sqrt(SNR_ul) * G * Xp + Wp;
% Despreading pilot 
Yp_aksen = Yp * fi;

% Channel Estimation -------------------------------------
% MMSE channel estimation
estG = zeros(M,K);            
for i=1:M;
    for j=1:K;
        % Estimated the chanel
        estG(i,j) = (sqrt(tau_p*SNR_ul)*beta) / (1+tau_p*SNR_ul*beta) .* Yp_aksen(i,j); 
    end
end
eEst = G - estG;              % Channel estimation error

%Mean-square channel estimation (unused, for verification only)
gamma = (tau_p*SNR_ul*beta^2) / (1+tau_p*SNR_ul*beta); 
MSE = beta - gamma;           % Mean Square Error

% Generate sample data for each UE -----------------------
jBit = nBit*K*2;              % Total generated bits
% Pre-allocating
b = zeros(K, jBit);
for i=1:K;
    % Generate streams of binary data (serial)
    b(i,:) = randi([0 1], [1 jBit]);     
end

% 4-QAM Modulation ---------------------------------------
symbol = QAM_symbol / sqrt(2); 
% Pre-allocating
q = zeros(K, jBit/nBit);
for i=1:K;
    % Divide streams into chunks (Serial to Parallel)
    pBit = transpose(reshape(b(i,:), nBit, []));
    % Convert binary into decimal
    pDec = bi2de(pBit);
    % Use decimal as 4 QAM inputs
    pComp = symbol(pDec+1, :);
    % Combine real and imaginer
    pQAM_real = transpose(pComp(:, 1));
    pQAM_imag = transpose(pComp(:, 2));
    pQAM = complex(pQAM_real, pQAM_imag);
    % QAM symbol for each UE
    q(i,:) = pQAM;
end

% Precoding ----------------------------------------------
% Identical power for each UE
power_dl = [1/K 1/K 1/K];     
% Zero Forcing Precoding Matrix
A = sqrt(M-K) * conj(estG) / (transpose(estG)*conj(estG));
% Signal to be transmitted (mixed)
x = A * sqrt(diag(power_dl)) * q; 

% OFDM ---------------------------------------------------
N = size(x, 1) ;              % Number of subcarrier = number of Tx antenna 
NFFT = size(x, 2)             % Number point of FFT must same as number of subcarrier
CP = ceil(0.25*NFFT);         % Number of cyclic Prefix (25% of NFFT)
% IFFT process
IFFT_SC = zeros(N, NFFT);
xCP = zeros(N, CP);
for i = 1:N;
    % IFFT of transmitted signal
    IFFT_SC(i,:) = ifft(x(i,:), NFFT);     
    % Copy the end of signal to the begining of signal
    xCP(i,:) = IFFT_SC(i, NFFT-CP+1:NFFT); 
end

% Add CP -------------------------------------------------
OFDM_symbol = [xCP IFFT_SC];  %CP + IFFT

% ===================== BS > UE ==========================
% Generate noise (AWGN)
W_real = rand(K, size(OFDM_symbol, 2));
W_imag = rand(K, size(OFDM_symbol, 2));
W = sqrt(vn) .* complex(W_real, W_imag) + mu;

% ===================== UE SIDE ==========================
% Received Signal at UE
HA = transpose(estG);
yk = sqrt(SNR_dl) * HA * OFDM_symbol + W;
 % Remove channel effect using zero forcing
yki = pinv(HA) * yk;

% Remove CP ----------------------------------------------
% Received OFDM symbol with noise and channel
OFDM_rem = yki(:, size(xCP, 2)+1:end);
fft_symbol = zeros(N, NFFT);
for i=1:N
    fft_symbol(i,:) = fft(OFDM_rem(i,:), NFFT);
end
xx = fft_symbol;

% Detection -----------------------------------------------
DA = transpose(estG);
Dq = DA*xx;

% 4-QAM De-Modulation ------------------------------------
yQAM = Dq;                     % Demodulation input
% Pre-allocating
yNumSymbol = size(yQAM, 2);
yb = zeros(K, nBit * yNumSymbol);
yDec = zeros(yNumSymbol, 1);
% Iterate every UE
for i = 1:K;
    % Separate real & imag from complex
    z9 = [real(yQAM(i,:)); imag(yQAM(i,:))];
    % Convert serial to parallel
    z0 = transpose(z9);
    % Get only the sign
    z1 = sign(reshape(z0, [], 2));
    % Iterate every symbol
    for j = 1:size(z1, 1);
        % Set true for matched QAM symbol
        z2 = ismember(QAM_symbol, z1(j, :), 'rows');
        % Get index of true result
        yDec(j) = find(z2) - 1;
    end;
    % Convert decimal to binary
    yBit = de2bi(yDec);  
    % Paralel to serial
    yb(i,:) = reshape(yBit', 1, []); 
end
BER = yb-b





