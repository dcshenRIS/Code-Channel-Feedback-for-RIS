function Rs=Rs_cal(H,H_RVQ,p)
%% descriptions of this function
% This function is to calculate the rate.
% ---------------- input descriptions -------------------------------------
%   "H" is the perfect channel. 
%   "H_RVQ" is the fed back channel.
%   "p" is the per-user transmit power
% ---------------- output descriptions ------------------------------------
%   "Rs" is the sum-rate.
M=size(H,1);
K=size(H,2);
W = H_RVQ/(H_RVQ'*H_RVQ);
W = W./repmat(sqrt(sum(abs(W).^2,1)),M,1);
Rs=0;
 for i_K = 1:K
     SINR = p*norm(H(:,i_K)'*W(:,i_K))^2/(1 + p*H(:,i_K)'*W*(H(:,i_K)'*W)'-p*norm(H(:,i_K)'*W(:,i_K))^2);   %p=gamma/K  
     Rs = Rs+log2(1+SINR);   
 end