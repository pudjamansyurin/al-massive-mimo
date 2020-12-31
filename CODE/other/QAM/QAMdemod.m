
symbolQAM = [-1 1; 1 1; 1 -1 ;-1 -1];
b = [3.6074    2.7006   -3.1977    3.9280
   -3.0065   -2.3885    1.2043   -4.1374
];
% Convert serial to parallel
bParallel = transpose(b);
% Get only the sign
bSign = sign(reshape(bParallel, [], 2));
% Iterate every symbol
bDecimal = zeros(1, size(bSign, 1));
for j = 1:size(bSign, 1);
    % Set true for matched QAM symbol
    sMatched = ismember(symbolQAM, bSign(j, :), 'rows')
    % Get index of true result
    bDecimal(j) = find(sMatched) - 1
end;
% Convert decimal to binary
bBinary = de2bi(bDecimal);  
% Paralel to serial
bSerial = reshape(transpose(bBinary), 1, []); 

