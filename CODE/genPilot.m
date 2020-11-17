function [ Xpf ] = genPilot( tau_p, K, N )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Generate a random complex matrix 
X = complex(rand(tau_p), rand(tau_p)) / sqrt(2);
% Factorize the matrix
[Q, R] = qr(X);
% Unitary matrix M x N
fi = Q(:, 1:K);
% Verification of unitary matrix
% verify = ctranspose(fi) * fi;
% Generate the pilot
Xp = sqrt(tau_p) * ctranspose(fi);
Xpf = fft(Xp,N,3);

end

