function [QuantH,err_CRVQ]= quantiz_CRVQ_real(H,R_sq,B)
%% descriptions of this function
% This function is to acquire the quan. channel accoring to the Statistics-based codebook feedback.
% ---------------- input descriptions -------------------------------------
%   "H" is the perfect channel with M(BS antennas) and K(users). 
%   "R_sq" is the square root of the  user’s channel correlation matrix.
%   "B" is the codebook size
% ---------------- output descriptions ------------------------------------
%   "QuantH" is is the fed back channel with codebook.

if B<=25
M = size(H, 1);
K = size(H, 2);
QuantH = zeros(M,K);
% Codebook
for i_K=1:1:K
C20 = randn(M,2^B) + 1j*randn(M,2^B);
C2 = R_sq*C20;
C2 = C2./repmat(sqrt(sum(abs(C2).^2,1)),M,1);  
h_k = H(:,i_K)/norm(H(:,i_K));
[temp1,index1]=max(real(h_k'*C2));
h_CRVQ = C2(:,index1)*norm(H(:,i_K));
QuantH(:,i_K)=h_CRVQ;
end
err_CRVQ = temp1;
else  fprintf('error');
err_CRVQ = 0;
end

