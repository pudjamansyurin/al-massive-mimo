clc;
clear all;
Nbps = 4;
NgType = 1 %CP
M = 4; %4QAM
Nfft=64; % FFT size
Ng=Nfft/4; % GI (Guard Interval) length (Ng=0 for no GI)
Nsym=Nfft+Ng; % Symbol duration
Nvc=Nfft/4; % Nvc=0: no VC (virtual carrier)
Nused=Nfft-Nvc;
Nframe=3; % Number of symbols per frame
X= randint(1,Nused*Nframe,M); % bit: integer vector
Xmod= qammod(X,M,0,'gray')/norm(Nbps);
x_GI= zeros(1,Nframe*Nsym) %Extende with CP
kk1=[1:Nused/2]; kk2=[Nused/2+1:Nused]; kk3=1:Nfft; kk4=1:Nsym;
for k=1:Nframe
X_shift= [Xmod(kk2) Xmod(kk1)]
x= ifft(X_shift);
x_GI(kk4)= guard_interval(Ng,Nfft,NgType,x);
kk1=kk1+Nused; kk2= kk2+Nused; kk3=kk3+Nfft; kk4=kk4+Nsym;
end

PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile ’dB’
Delay=[0 3 5 6 8]; % Channel delay ’sample’
Power=10.^(PowerdB/10); % Channel tap power profile ’linear scale’
Ntap=length(PowerdB); % Chanel tap number
channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);
Lch=Delay(end)+1; % Channel length
h=zeros(1,Lch); h(Delay+1)=channel; % cir: channel impulse response
y = conv(x_GI,h);
kk1=(NgType==2)*Ng+[1:Nsym]; kk2=1:Nfft;
kk3=1:Nused; kk4=Nused/2+Nvc+1:Nfft; kk5=(Nvc~=0)+[1:Nused/2];
H= fft([h zeros(1,Nfft-Lch)]);
H_shift(kk3)= [H(kk4) H(kk5)];



% Coba sendiri
a = [1; 2; 3; 4; 5; 6; 0; 0;]
iffta = ifft(a)
CP = iffta(1:2)
ifftcp = vertcat(iffta,CP)
ffta = fft(ifftcp)
