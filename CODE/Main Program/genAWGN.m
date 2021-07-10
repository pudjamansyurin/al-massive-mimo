function [ nt, nf ] = genAWGN ( X, Y, N )

nt = sqrt(0.5)*(randn(X,Y,N)+1i*randn(X,Y,N)); % time domain
nf = sqrt(1/N)*fft(nt,N,3); % frequency domain

end

