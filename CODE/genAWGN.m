function [ nt, nf ] = genAWGN ( X, Y, N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nt = sqrt(0.5)*(randn(X,Y,N)+1i*randn(X,Y,N)); % time domain
nf = sqrt(1/N)*fft(nt,N,3); % frequency domain

end

