clc;
clear all;

K = 16;
M = [20:5:100];
SNR_DL = 1; %0dB
Beta = 1;
power = 1/K;
f = 30*10^9;
Bc = 20*10^6;
v = 142;
c = 3*10^8;

lambda = c/f;
Tc = lambda/(2*v);
tau_c = Bc*Tc;
tau_p = 0.2*tau_c; %asumsi pilot 20% 

%Zero Forcing;
for i = 1:length(M);
    C(i) = 16*(log(1+(M(i)-K)*SNR_DL*Beta*power))/log(2); %Zero forcing capacityZero forcing capacity
    Spectral(i) = (1-(tau_p/tau_c))*C(i);
end
%Maximum Ratio
for i = 1:length(M);
    CM(i) = 16*(log(1+(M(i)*SNR_DL*Beta*power)/(1+SNR_DL*Beta*16*power))/log(2));
    SpectralM(i) = (1-(tau_p/tau_c))*CM(i);
end
%zero forcing
figure(1)
plot(M,C) 
hold on
plot(M,CM)
title('Kapasitas Total Sistem Massive MIMO ');
xlabel('Jumlah antena BTS (M)');
ylabel('Kapasitas');
legend('Zero Forcing','Maximum Ratio');
grid on;
%Maximum Ratio
figure(2)
plot(M,Spectral)%zero forcing spectral efficiency
hold on;
plot(M,SpectralM);
title('Spectral Efficiency Total Sistem Massive MIMO');
xlabel('Jumlah antena BTS (M)');
ylabel('Spectral Effifiency ');
legend('Zero Forcing','Maximum Ratio');
grid on;
figure(3)
plot(M,C/16); 
hold on
plot(M,CM/16);
title('Kapasitas Tiap User Sistem  Massive MIMO');
xlabel('Jumlah antena BTS (M)');
ylabel('Kapasitas');
legend('Zero Forcing','Maximum Ratio');
grid on;
figure(4)
plot(M,Spectral/16)%zero forcing spectral efficiency
hold on;
plot(M,SpectralM/16);
title('Spectral Efficiency Tiap User Sistem Massive MIMO');
xlabel('Jumlah antena BTS (M)');
ylabel('Spectral Effifiency ');
legend('Zero Forcing','Maximum Ratio');
grid on;







