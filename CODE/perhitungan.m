close all;
clc;
clear all;
teta = [-30 10 40];
% teta = [-80 -65 -45 -30 -10 5 15 40 55 85 ]
M = 300;
gamma = -90:0.1:90;
dH = 1/2;
AF = zeros(length(teta),length(gamma));
normAF = zeros(length(teta),length(gamma));
for tt = 1:length(teta);
    num = pi*dH*M*(sind(teta(tt))-sind(gamma));
    denum = pi*dH*(sind(teta(tt))-sind(gamma));
    AF(tt,:) = sind(num).^2./(M.*(sind(denum)).^2);
    normAF(tt,:) = AF(tt,:)./max(AF(tt,:));
    plot(gamma,normAF(tt,:));
%     polar(gamma*pi/180, normAF(tt,:));
    hold on;
end
% plot(gamma,sum(normAF,1))
% polar(gamma*pi/180, sum(normAF,1));

% legend(sprintfc('Sudut UE = %d',teta));
title('Array Response Massive MIMO (M = 300, \theta = 30^{0})');
% title('Array Response Massive MIMO (M = 300)');
% legend(fprintf('\theta = %d', teta))
% legend(sprintfc('\theta ', teta));
xlabel('Sudut UE (\theta)');
ylabel('Normaliasi Array Response');
% title('Array Respons Kondisi LOS');


% close all;
% clc;
% clear all;
% teta =[-40 10];
% M = 300;
% gamma = -180:1:180;
% dH = 1/2;
% AF = zeros(length(teta),length(gamma));
% normAF = zeros(length(teta),length(gamma));
% % for tt = 1:length(teta);
%     num = pi*dH*M*(sind(teta(1))-sind(teta(2))-sind(gamma));
% %     denum = pi*dH*(sind(teta(tt))-sind(teta(tt+1))-sind(gamma));
% %     AF = sind(num).^2./(M.*(sind(denum)).^2);
% %     normAF(tt,:) = AF(tt,:)./max(AF(tt,:));
% % %     plot(gamma,normAF(tt,:));
% % %     polar(gamma*pi/180, normAF(tt,:));
% % % %     hold on;
% % end









