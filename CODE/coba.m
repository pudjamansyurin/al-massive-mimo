% M = 20:10:100;
% K = 10;
% teta = [10 20 30];
% gamma = teta+0.001;
% dH = 1/2;
% for m = 1:length(M);
%     for t = 1:length(teta)
%         num = pi*dH*M(m)*(sind(teta(t))-sind(gamma));
%         denum = pi*dH*(sind(teta(t))-sind(gamma));
%         AF(t,:) = sind(num).^2./(M(m).*(sind(denum)).^2);
%     end
%     sumAF = mean(sum(AF,2)-diag(AF))./max(AF);
% end

Ko = 6;
 for k = 1:Ko
     sinteta(k) = -1+(2*k-1)/Ko;
 end