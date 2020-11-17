par = getParam();

% validation
if length(par.K) > 1
    fprintf('K => should scalar!\n');
    return
end

% generate sample data
xb = genData(par.FRM, par.BPU);

rBER = zeros(length(par.Code), length(par.SNR_dB)); 
for Si = 1 : length(par.SNR_dB); 
    SNRo = par.SNR_dB(Si);
    N0 = 10.^(-SNRo/10);

    % Loop through all data ----------------------------------
    yb = zeros(length(par.Code), par.FRM, par.BPU);
    for Fi=1:par.FRM;

        % 4-QAM Modulation ------------------------------------
        q = modQAM(par.BPS, par.symbol, xb(Fi,:));

        % Convert 2D matrix into 3D matrix
%         S = reshape(q, size(q,1), 1, size(q,2));

        % Channel ---------------------------------------
        [Ht, Hf] = genChannel(Cho,Ko,L,Mo,N,beta,teta);

        % Channel estimation ---------------------------
        [ nt, nf ] = genAWGN ( Mo, tau_p, N );                        
        Xpf = genPilot( tau_p, K, N );
        Hf_est = zeros(Ko,Mo,N);
        Error = zeros(Ko,Mo,N);

        for n = 1:N;  
            Hf_est(:,:,n) = estChannel( Hf(:,:,n), Xpf(:,:,n), N0, nf(:,:,n) );
            Error(:,:,n) = Hf_est(:,:,n) - Hf(:,:,n);
            MSE(Si) = mean(mean(abs(Error(:,:,n).^2)));
        end

%                     Hf = Hf_est;

        % Generate AWGN noise
        [nt, nf] = genAWGN(Ko, 1, N);

        % Loop multiple precode --------------------------------
        meanRk = zeros(length(Code),Ko,1);
        for Ci=1:length(Code);

            % Precode
            A = zeros(Mo,Ko,N);
            C = zeros(Mo,1,N);
            for i = 1:N;
                [A(:,:,i), C(:,:,i)] = genPrecoding(Code(Ci), Ko, Mo, Hf(:,:,i), S(:,:,i),N,Pc);
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
                [SINR(:,:,n),Sig,I,Noise,Pt,AA,Sig2] = calcSINR(A(:,:,n), Hf(:,:,n), SNRo, Ko, Pc);
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

        end % End of precoding loop

    end % End of frames loop

    % BER calculation
    xbr = reshape(xb, [1, numel(B)]);
    for Ci=1:length(Code);
        ybr = reshape(yb(Ci,:,:,:), [1, numel(B)]);
        [eBER(Ci,Si), rBER(Ci,Si)] = biterr(xbr, ybr);
    end
end % End of SNR loop

% Plot Bit Error Rate MU-Massive MIMO
figure(1);
for Ci=1:length(par.Code);
    semilogy(par.SNR_dB, rBER(Ci,:), genMark(1,Ci,1));
    hold on;
end
grid on;
legend(strcat({''}, par.Code));
title('Bit Error Rate MU-Massive MIMO');
xlabel('SNR(dB)');
ylabel ('BER');