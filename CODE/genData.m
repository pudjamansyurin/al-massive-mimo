function [ xb ] = genData( frame, bit_per_user )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Generate sample data for each UE -----------------------
xb = randi([0, 1], [frame, bit_per_user]);
            
end

