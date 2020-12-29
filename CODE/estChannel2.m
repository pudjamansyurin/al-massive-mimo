function [ Hf_est,Yp ] = estChannel2( Xp,Ht,M,K,tau_p )

% nois = sqrt(0.5)*(randn(M,tau_p)+1i*randn(M,tau_p))
Yp = Ht'*Xp;
Hf_est = Yp/Xp;


% estG = zeros(M,K);            
% for i=1:M;
%     for j=1:K;
%         % Estimated the chanel
%         Hf_est(i,j) = (sqrt(tau_p*SNR_L)*beta) / (1+tau_p*SNR_L*beta) .* Yp(i,j); 
%     end
% end
%              % Channel estimation error


end

