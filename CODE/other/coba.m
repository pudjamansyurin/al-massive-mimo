clear all; clc;         % Clear screen
tau_p = 5;              % Pilot length
M = 10;                 % Number of Tx antenna (in one BS)
K = 3;                  % Number of Rx antenna (= number of UE)
L = 4;                  % Channel tap frequency selective
vn = 1;                 % Variance noise
mu = 0;                 % Mean of noise
beta = 0.5;               % Large scale fadding Coefficient Of Rayleigh Channel
% SNR_ul = 128;           % Uplink SNR (dB)
% SNR_dl = SNR_ul;        % Download SNR (dB)
BPS = 2;                % (Bit/Symbol) Number of bits per chunk
N = 16;                 % Number of subcarrier
NFFT = L*N;             % Number point FFT  
nBit = 2;
QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];

% Generate sample data for each UE -----------------------
BPU = N*2;              % (Bit/User)    
% Pre-allocating
b = zeros(K, BPU);
for i=1:K;
    % Generate streams of binary data (serial)
    b(i,:) = randi([0 1], [1 BPU]);     
end

% 4-QAM Modulation ---------------------------------------
symbol = QAM_symbol / sqrt(2); 
% Pre-allocating
q = zeros(K, BPU/BPS);
for i=1:K;
    % Serial to Parallel (by BitPerSymbol)
    pBit = transpose(reshape(b(i,:), BPS, []));
    % Convert binary into decimal (start from 1, not 0)
    pDecimal = bi2de(pBit) + 1;
    % Use decimal as 4 QAM inputs
    pComplex = symbol(pDecimal, :);
    % Combine real and imaginer
    % pQAM_real = transpose(pComplex(:, 1));
    % pQAM_imag = transpose(pComplex(:, 2));
    pQAM = complex(pComplex(:, 1)', pComplex(:, 2)');
    % QAM symbol for each UE
    q(i,:) = pQAM;
    
%    
end
sq = size(q);
S = reshape(q, sq(1), 1, sq(2));


% ================== Rayleigh Channel ====================
% Standard deviation
Pch_DB = [-10 -10 -10 -10];         % Channel tap power profile ’dB’
Pch = 10.^(Pch_DB/10);              % Channel tap power profile ’linear scale’
v = 1./sqrt(Pch);                   % Variance channel for frequency selective
sigma = sqrt(v);
% Generate the Rayleigh channel from UE to BS
H = zeros(M,K,L);
G = zeros(M,K,L);
% -------------------------New-------------------
% Time-domain Rayleigh channel matrix BTS side
  Ht = sqrt(0.5/L)*(randn(M,K,L) + 1i*randn(M,K,L));
  Hta = sqrt(0.5/L)*(randn(K,M,L) + 1i*randn(K,M,L));
  Hfa = fft(Hta,N,3);
% New Precoding
for i = 1:N;
    A(:,:,i) = Hfa(:,:,i)'/(Hfa(:,:,i)*Hfa(:,:,i)'); % precoding matrix
end
%-----------------------End New-----------------
% Frequency domain channel for each subcarrier

% QAM * Precoding ----------------------------------------
% Pre-allocating
C = zeros(M,1,N);
for i=1:N;
    C(:,:,i) = A(:,:,i)*S(:,:,i);
end

%----------RECEIVER---------------
nt = sqrt(0.5)*(randn(par.U,par.N,par.nrof_ofdm_symbols)+1i*randn(par.U,par.N,par.nrof_ofdm_symbols)); % time domain
 nf = sqrt(1/par.N)*fft(nt,par.N,2); % frequency domain
for i = 1:N;  
    % Received Signal
    yf(:,:,i) = Hfa(:,:,i)*C(:,:,i);   
end

At = yf(:,:,:);
xt = size(At);
Bt = reshape(At, xt(1), xt(3));

yQAM = Bt;
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