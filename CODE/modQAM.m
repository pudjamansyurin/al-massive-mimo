function [ QAM ] = modQAM( BPS, symbol, binary_data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Serial to Parallel (by BitPerSymbol)
parallel_bit = transpose(reshape(binary_data, BPS, []));
% Convert binary into decimal (start from 1, not 0)
decimal_data = bi2de(parallel_bit) + 1;
% Use decimal as 4 QAM inputs
complex_value = symbol(decimal_data, :);
% Combine real and imaginer
QAM = complex(complex_value(:, 1)', complex_value(:, 2)');

end

