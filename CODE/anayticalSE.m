Analytical SE

        Anlytical Spectral Efficiency each user (ZF)
        SINRt(Mi) = (Mo-Ko)*SNR_L*beta*Pc;    % SINR theory each user
        SEt(Mi) = log2(1+SINRt(Mi));
        SEt_all(Mi) = Ko*SEt(Mi);
%         SINRt(Mi) = Pc*beta^2*Ko^2/((Ko*beta-beta)*Pc*beta*Mo + (beta^2*Mo*Ko^2));
%         SINRt(Mi) = 10*(Mo-Ko)/Ko;
%         SEt_all(Mi) = Ko*log2(1+SINRt(Mi));
        
        % Anlytical Spectral Efficiency each user (MRT)
        SINRtm(Mi) = Mo*SNR_L*beta*Pc/(1+SNR_L*beta*Ko*Pc);
        SEtm(Mi) = log2(1+SINRtm(Mi));
        SEtm_all(Mi) = Ko*SEtm(Mi);
        
        % Analytical NLOS
        SINRp(Mi) = (Mo-1)/((Ko-1)+Ko*0.001+(1/SNR_L));
        SEtp(Mi) = Ko*log2(1+SINRp(Mi));
            
        % Analytical LOS
        SINR_LOS(Mi) = 1/(beta + 1/SNR_L);
        SE_LOS(Mi) = K*log2(1+SINR_LOS);

plot(M,SEt_all,'--m');
hold on;
plot(M,SEtm_all,'--c');
hold on;
plot(M,SEtp,'--g');
hold on;
plot(M,SE_LOS,'--g');