function [ xb ] = genData( frame, bit_per_user )

% Generate sample data for each UE -----------------------
xb = randi([0, 1], [frame, bit_per_user]);
            
end

