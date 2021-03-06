function [ Ht, Hf ] = genChannel(type, K, L, M, N,beta,teta)

%======================UR-LOS channel==================
if strcmp(type, 'LOS')
    comp = zeros(K,M);
    Ht = zeros(K,M);
    for k = 1:K
        for m = 1:M;
            comp(k,m) = exp(-1i*(m-1)*pi*sind(teta(k)));
        end
        % Time-domain LOS channel matrix BTS side
        Ht = sqrt(beta)*comp; 
    end
    % Freq-domain LOS channel matrix BTS side
    Hf = fft(Ht,N,3); 
%======================Rayleigh Channel================
else
    % Time-domain Rayleigh channel matrix BTS side
    Ht = sqrt(0.5/L)*(randn(K,M,L) + 1i*randn(K,M,L));
%     scale=1e-6; % nano
%     sigma_t=2.3*scale; % RMS delay spread
%     Ts = 10*sigma_t/L;
%     sigma02=(1-exp(-Ts/sigma_t))/(1-exp(-(L+1)*Ts/sigma_t)); % Eq.(2.9) 
%     l=0:L-1; 
%     PDPo = sigma02*exp(-l*Ts/sigma_t);
%     for l = 1:L;
%         Ht(:,:,l) = sqrt(PDPo(l)/2)*((randn(K,M) + 1i*randn(K,M)));
%     end
    % Freq-domain Rayleigh channel matrix BTS side
    Hf = fft(Ht,N,3);
end