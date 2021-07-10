function [bSerial ] = demodQAM( b, symbolQAM )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Convert serial to parallel
bParallel = transpose(b);
% Get only the sign
bSign = sign(reshape(bParallel, [], 2));
% Iterate every symbol
bDecimal = zeros(1, size(bSign, 1));
for j = 1:size(bSign, 1);
    % Set true for matched QAM symbol
    sMatched = ismember(symbolQAM, bSign(j, :), 'rows');
    % Get index of true result
    bDecimal(j) = find(sMatched) - 1;
end;
% Convert decimal to binary
bBinary = de2bi(bDecimal);  
% Paralel to serial
bSerial = reshape(transpose(bBinary), 1, []); 

end

