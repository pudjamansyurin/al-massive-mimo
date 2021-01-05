% Clear screen
clear all; close all; clc;         

% M = 50;               % Number of Tx antenna (in one BS)
M = 50:50:350;              % Number of occupied subcarrier
N = 350;                 % Number of Rx antenna (= number of UE)
K = 30;
L = 4;                  % Channel tap frequency selective
beta = 1;
BPS = 2;                % (Bit/Symbol) Number of bits 
nBit = 2;               % Numer \bit per symbol
nCP = ceil(0.05*N);     % Number of cyclic Prefix (25% of NFFT)
% SNR_dB = -10:2:20;    % list of SNR [dB] values to be simulated
SNR_dB = 0;
SNR_L = 10^(SNR_dB(length(SNR_dB))/10);
FRM = 1;              % Number of data frame
tau_p = 10;
BPU = N*2;              % (Bit/User)  
NBPU = BPU*FRM;

QAM_symbol = [-1 1; 1 1; 1 -1 ;-1 -1];
symbol = QAM_symbol / sqrt(2); 
Code = { 
     'ZF' 
     'MRT' 
     'MMSE' 
};
Channel = { 
    'LOS' 
%     'Rayleigh' 
};
CSI = {
    'Perfect CSI'
%     'Imperfect CSI'    
};

%  ===================START SIMULATION====================
fprintf('running numerical simulation. \n');
%====================TRANSMITTER SIDE=====================
Prog = 0;
OneProg = (100/(length(CSI)*length(K)*length(M)*length(Channel)*length(SNR_dB)*FRM*length(Code)));

% Loop CSI --------------------------
SE = zeros(length(CSI), length(Channel), length(K), length(Code), length(M));
rBER = zeros(length(CSI), length(Channel), length(Code), length(SNR_dB)); 
for Ei = 1:length(CSI);
    Eo = CSI(Ei);
    
    % Loop rx antenna --------------------------
    for Ki = 1:length(K);
        Ko = K(Ki);
        Pc = 1/K(Ki);   % Power control each user
        for k = 1:Ko
            sinteta(k) = -1+(2*k-1)/Ko;
        end
        teta = asind(sinteta);
        
        % Loop tx antenna --------------------------
        SE_LOS = zeros(1,length(M));
        SE_NLOS = zeros(1,length(M));
        SE_NZF = zeros(1,length(M));
        for Mi = 1:length(M);
            Mo = M(Mi);

            % Loop Channel--------------------------
            MSE = zeros(length(Channel), length(SNR_dB));
            for Chi = 1:length(Channel)
                Cho = Channel(Chi);

                % Generate sample data for each UE -----------------------
                B = randi([0, 1], [Ko NBPU]);
                % Divide into frames
                xb = reshape(B, [FRM, Ko, BPU]);

                % Loop through all SNR ----------------------------------
                for Si = 1 : length(SNR_dB); 
                    SNRo = SNR_dB(Si);
                    N0 = 10.^(-SNRo/10);

                    % Loop through all data ----------------------------------
                    yb = zeros(length(Code), FRM, Ko, BPU);
                    for Fi=1:FRM;

                        % 4-QAM Modulation ---------------------------------------
                        q = zeros(Ko, BPU/BPS);
                        for k=1:Ko;
                            % QAM symbol for each UE
                            q(k,:) = modQAM(BPS, symbol, xb(Fi,k,:));
                        end

                        % Convert 2D matrix into 3D matrix
                        S = reshape(q, size(q,1), 1, size(q,2));

                        % Channel ---------------------------------------
                        [Ht, Hf] = genChannel(Cho,Ko,L,Mo,N,beta,teta);

                        % Channel estimation ---------------------------
                        if strcmp(Eo, 'Imperfect CSI')
                            [ nt, nf ] = genAWGN ( Mo, tau_p, N );                        
                            [Xpf,Xp,fi] = genPilot( tau_p, K, N );
                            
                            Hf_est = zeros(Ko,Mo,N);
                            for n = 1:N;  
                                fi2 = fft(fi,N,3);
                                Hf_est(:,:,n) = estChannel( Hf(:,:,n), Xpf(:,:,n), N0, nf(:,:,n), tau_p, SNRo,fi2(:,:,n));
                            end
                            
                            Error = abs((Hf_est - Hf).^2);
                            MSE(Chi,Si) = mean(reshape(Error, [1 numel(Error)]));
                            
                            Hf = Hf_est;
                        end;

                        % Generate AWGN noise
                        [nt, nf] = genAWGN(Ko, 1, N);

                        % Loop multiple precode --------------------------------
                        meanRk = zeros(length(Code),Ko,1);
                        for Ci=1:length(Code);

                            % Precode
                            A = zeros(Mo,Ko,N);
                            C = zeros(Mo,1,N);
                            for i = 1:N;
                                [A(:,:,i), C(:,:,i)] = genPrecoding(Code(Ci), Ko, Mo, Hf(:,:,i), S(:,:,i),N,Pc,SNRo);
                            end

                            % Transform to time domain
                            xt = ifft(C,N,3); % transform to time domain

                            % Cyclic Prefix
                            CP = zeros(Mo,1,nCP);
                            for i = 1:Mo;
                                CP(i,1,:) = xt(i, :, N-nCP+1:length(xt));
                            end

                            % Combine CP with signal
                            xt_CP = cat(3, CP, xt);

                            %===================RECEIVER SIDE=========================
                            % Remove cyclic prefix
                            xt_rem = xt_CP(:,1, nCP+1:length(xt_CP));

                            % Transform to freq domain
                            xf = fft(xt_rem,N,3);                        

                            % Multipliying channel and signal for each subcarrier
                            Hfy = zeros(Ko,1,N);
                            for i = 1:N;  
                                Hfy(:,:,i) = Hf(:,:,i)*xf(:,:,i);
                            end

                            % noise variance
                            yf = Hfy + sqrt(N0)*nf;

                            % Detection (Demodulasi QAM)
                            yQAM = reshape(yf, size(yf,1), size(yf,3));
                            for i = 1:Ko;
                                % Separate real & imag from complex
                                qm = [real(yQAM(i,:)); imag(yQAM(i,:))];
                                % Demodulating
                                yb(Ci,Fi,i,:) = demodQAM(qm, QAM_symbol);
                            end

                            % SINR--------------------------------------------------------
                            SINR = zeros(Ko,1,N);
                            sigmaA = zeros(Mo,Mo,N);
                            for n = 1:N;
                                [SINR(:,:,n),Sig,I,Noise,Pt,AA,Sig2] = calcSINR(A(:,:,n), Hf(:,:,n), SNR_dB, Ko, Pc);
                            end

                            %reshape 3D to 2D matrix
                            SINRa = reshape(SINR, size(SINR,1), size(SINR,3));

                            %Spectral Efficiency for each subcarrier
                            Rk = zeros(Ko,N);

                            for k = 1:Ko;
                                for n = 1:N;
                                    %Spectral Efficiency at all user and all subcarrier
                                    Rk(k,n) = log2(1 + SINRa(k,n)); 
                                end
                                % Mean Spectral Efficiency of all subcarrier
                                meanRk(Ci,k,:) = mean(Rk(k,:));
                            end

                            % Display info
                            Prog = Prog + OneProg;
                            fprintf('RUN: %.1f %%\n', Prog);
                        end % End of precoding loop

                    end % End of frames loop

                    % BER calculation
                    xbr = reshape(xb, [1, numel(B)]);
                    for Ci=1:length(Code);
                        ybr = reshape(yb(Ci,:,:,:), [1, numel(B)]);
                        [eBER, rBER(Ei,Chi,Ci,Si)] = biterr(xbr, ybr);
                    end
                end % End of SNR loop

                % Sum spectral Efficiency all user
                for Ci=1:length(Code);
                    SE(Ei, Chi, Ki,Ci,Mi) = sum(meanRk(Ci,:,:));
                end
            end  
             % =======================Analytical SE LOS=========================
        dH = 1/2;
        gamma = teta+0.001;
        for tt = 1:length(teta);
            for gm = 1:length(gamma);
                num = pi*dH*Mo*(sin(gamma(gm))-sin(teta(tt)));
                denum = pi*dH*(sin(gamma(gm))-sin(teta(tt)));
                AF(tt,gm) = sin(num).^2./(Mo*sin(denum).^2);        
            end
            AFsum = (sum(AF,2)-diag(AF));
        end
        AFsum(Mi) = mean(sum(AF,2)-diag(AF)); 
        AFsum(Mi) = AFsum(Mi)/max(AFsum(Mi));
        SE_LOS(Mi) = Ko*log2(1+(Mo./(AFsum(Mi)+1/SNR_L)));
         % =======================Analytical SE NLOS=========================
        SE_NLOS(Mi) = Ko*log2(1+((Mo-1)./((Ko-1)*(Mo-1)/Mo+Ko+1/SNR_L)));
        % =======================Analytical SE NLOS zf======================
        SE_NZF(Mi) = Ko*log2(1+Pc*(Mo-Ko)*SNR_L);
         end % End of tx antenna loop
    end % End of user loop
end

% Plot Bit Error Rate MU-Massive MIMO
figure(1);
lgd = cell(length(Channel)*length(Code),1);
i=1;
for Ei=1:length(CSI);
    for Chi=1:length(Channel);
        for Ci=1:length(Code);
            REE = rBER(Ei,Chi,Ci,:);
            RER = reshape(REE, [1, numel(REE)]);     
            semilogy(SNR_dB, RER, genMark(Chi,Ci,Ei));
            hold on;
            lgd(i) = strcat(CSI(Ei), {', '}, Channel(Chi), {', '}, Code(Ci)); 
            i = i + 1;
        end
    end 
end
grid on;
legend(lgd);
% title('Bit Error Rate (BER) MU-Massive MIMO');
xlabel('SNR(dB)');
ylabel ('Bit Error Rate (BER)');
        
% Plot Spectral Efficiency MU-Massive MIMO (Precoding)
figure(2);
% lgd = cell(length(Channel)*length(Code),1);
i = 1;
for Ei=1:length(CSI);
    for Chi=1:length(Channel);
        for Ci=1:length(Code);
            SEE = SE(Ei,Chi,Ki,Ci,:);
            SEX = reshape(SEE, [1, numel(SEE)]);
            plot(M, SEX, genMark(Chi, Ci, Ei));
            hold on;    
            lgd(i) = strcat(CSI(Ei), {', '}, Channel(Chi), {', '}, Code(Ci)); 
            i = i + 1;
        end    
    end
end
plot(M,SE_NLOS,'--')
hold on;
% plot(M,SE_NZF,'-*')
lgd(i) = strcat({'Lower bound NLOS'}); 
legend(lgd);
grid on;
% title(sprintfc('Spectral Efficiency MU-Massive MIMO (K = %d)', K));
xlabel('Jumlah antena BTS (M)');
ylabel('Efisiensi spektrum (Bit/s/Hz)');
    
% Plot Spectral Efficiency MU-Massive MIMO (User)
figure(3);
Coding = find(ismember(Code, 'MMSE'));
for Ki = 1:length(K);
    SEE = SE(1,Chi,Ki,Coding,:);
    SEX3 = reshape(SEE, [1, numel(SEE)]);
    plot(M, SEX3, genMark(1,Ki,1));
    hold on;
end
grid on;
legend(sprintfc('K = %d', K));
% title('Spectral Efficiency MU-Massive MIMO');
xlabel('Jumlah antena BTS (M)');
ylabel('SEfisiensi spektrum (Bit/s/Hz)');

% Plot Spectral Efficiency MU-Massive MIMO (User)
figure(4);
Coding = find(ismember(Code, 'MMSE'));
for Mi = 1:length(M);
    SEE = SE(1,Chi,:,Coding,Mi);
    SEX4 = reshape(SEE, [1, numel(SEE)]);
    plot(K, SEX4, genMark(1,Mi,1));
    hold on;
end
grid on;
legend(sprintfc('M = %d', M));
% title('Spectral Efficiency MU-Massive MIMO');
xlabel('Jumlah User (K)');
ylabel('Efisiensi spektrum (Bit/s/Hz)');

figure(5);
for Chi=1:length(Channel);
    plot(SNR_dB, MSE(Chi,:), genMark(Chi,Chi,1));
    hold on;
end
grid on;
legend(strcat({''}, Channel));
% title('MSE of Channel Estimation');
xlabel('SNR (dB)');
ylabel('Mean-Squared Error (MSE)');

figure(6);
plot(M,SE_LOS);
hold on;
plot(M,SE_NLOS);
legend('LOS','NLOS');
title('Analytical Spectral Efficiency');
xlabel('Jumlah Antena BTS (M)');
ylabel('Spectral Efficiency (Bit/s/Hz)');

% % Save to log
% filename = strcat('logs/', datestr(now, 30), '.mat');
% save(filename);


    