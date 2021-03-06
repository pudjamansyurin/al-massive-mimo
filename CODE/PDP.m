scale=1e-6 % nano
sigma_t=4.7*scale; % RMS delay spread
L = 22;
Ts = 10*sigma_t/L;
% Ts = 50*scale;
% L = 10*sigma_t/Ts
sigma02=(1-exp(-Ts/sigma_t))/(1-exp(-(L+1)*Ts/sigma_t)); % Eq.(2.9)
l=0:L; 
PDPo = sigma02*exp(-l*Ts/sigma_t); % Eq.(2.8)
stem(l,PDPo);
xlabel('Indeks Delay Tap Kanal')
ylabel('Power Delay Profile (PDP)')