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
    % Freq-domain Rayleigh channel matrix BTS side
    Hf = fft(Ht,N,3);
end