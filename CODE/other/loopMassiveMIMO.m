clear all; clc;         % Clear screen
% M = 10;                 % Number of Tx antenna (in one BS)
%num = 100:100:500;     % Number antenna Tx
K = 3;                 % Number of Rx antenna (= number of UE)
L = 4;                  % Channel tap frequency selective
beta = 0.5;             % Large scale fadding Coefficient Of Rayleigh Channel
BPS = 2;                % (Bit/Symbol) Number of bits per chunk
N = 16;                % Number of occupied subcarrier
num = 20:10:30;
nBit = 2;               % Number bit per symbol
QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];
SNR_dB = 0:1:20;        % list of SNR [dB] values to be simulated

% Loop antenna Tx----------------------------------------
ber.error = zeros(length(num), length(SNR_dB));
ber.ratio = zeros(length(num), length(SNR_dB)); 

for tx = 1:length(num);
    progressTX = (tx-1)*100/length(num);
%     L = num(tx);
    M = num(tx);
    nCP = ceil(0.25*N);     % Number of cyclic Prefix (25% of NFFT)
    %====================TRANSMITTER SIDE=====================
    % Generate sample data for each UE -----------------------
    BPU = N*2;              % (Bit/User)  
    NBPU = BPU*10;
    chunk = NBPU/BPU;

    B = randi([0, 1], [K NBPU]);
    % Divide into chunks
    xb = reshape(B, [chunk, K, BPU]);

    % Loop through all SNR ----------------------------------
    for snr = 1 : length(SNR_dB); 
        progressBER = (snr-1)*100/length(SNR_dB);
            
        % Loop through all data ----------------------------------
        yb = zeros(chunk, K, BPU);
        for X=1:chunk;
            progressChunk = (X*BPU)*100/NBPU;

            % 4-QAM Modulation ---------------------------------------
            symbol = QAM_symbol / sqrt(2); 
            q = zeros(K, BPU/BPS);
            for i=1:K;
                % QAM symbol for each UE
                q(i,:) = modQAM(BPS, symbol, xb(X,i,:));
            end
            % Convert 2D matrix into 3D matrix
            S = reshape(q, size(q,1), 1, size(q,2));

            % Rayleigh Channel ---------------------------------------
            Hf = generateChannel(K,L,M,N);

            [A, C] = generatePrecoding(K,M,N,Hf,S);

            % Transform to time domain
            xt = ifft(C,N,3); % transform to time domain

            % Cyclic Prefix
            CP = zeros(M,1,nCP);
            for i = 1:M;
                CP(i,1,:) = xt(i,:,N-nCP+1:length(xt));
            end
            % Convert 2D to 3D matrix
            xt_2D = reshape(xt, size(xt,1), size(xt,3));
            CP_2D = reshape(CP, size(CP,1), size(CP,3));

            % Add cyclic perfix
            xt_CP_2D = [CP_2D xt_2D];

            % Convert 3D to 2D matrix (Time domain signal with cyclic prefix)
            xt_CP = reshape(xt_CP_2D,size(xt_CP_2D,1), 1, size(xt_CP_2D,2)); 

            %===================RECEIVER SIDE=========================
            % Remove cyclic prefix
            xt_rem = xt_CP(:,1, nCP+1:length(xt_CP));
            % Transform to freq domain
            xf = fft(xt_rem,N,3);
            % Generate noise
            nt = sqrt(0.5)*(randn(K,1,N)+1i*randn(K,1,N)); % time domain
            nf = sqrt(1/N)*fft(nt,N,3); % frequency domain

            Hfy = zeros(K,1,N);
            for i = 1:N;  
                % Multipliying channel and signal for each subcarrier
                Hfy(:,:,i) = Hf(:,:,i)*xf(:,:,i);   
            end
            % noise variance    
            N0 = 10.^(-SNR_dB(snr)/10);
            yf = Hfy + sqrt(N0)*nf;
            % Reshape 3D matrix into 2D matrix
            At = yf(:,:,:);
            Bt = reshape(At, size(At,1), size(At,3));
            % Detection (Demodulasi QAM)
            yQAM = Bt;
            for i = 1:K;
                % Separate real & imag from complex
                qm = [real(yQAM(i,:)); imag(yQAM(i,:))];
                % Demodulating
                yb(X,i,:) = demodQAM(qm, QAM_symbol);
            end
            % Display info
%             progress = progressTX + progressBER/length(num) + (progressChunk/length(SNR_dB));
%             fprintf('RUN: %.3f %%\n', progress);
            fprintf('TX: %.1f %%, ', progressTX);
            fprintf('SNR: %.1f %%, ', progressBER);
            fprintf('CHUNK: %.1f %%\n', progressChunk);
        end % End of chunk loop

        % BER calculation
        xbr = reshape(xb, [1, numel(B)]);
        ybr = reshape(yb, [1, numel(B)]);
        [ber.error(tx, snr), ber.ratio(tx, snr)] = biterr(xbr, ybr);

    end % End of SNR loop
end

semilogy(SNR_dB, ber.ratio)
legend('N = 300','N = 400','N = 500');