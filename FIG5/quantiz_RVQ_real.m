function [QuantH,err_RVQ]= quantiz_RVQ_real(H, B)
% only quantize the row direction
if B<=25
M = size(H, 1);
K = size(H, 2);
QuantH=zeros(size(H));
C1=randn(M,2^B) + 1j*randn(M,2^B);
C1 = C1./repmat(sqrt(sum(abs(C1).^2,1)),M,1);  
for i_K=1:1:K
h_k = H(:,i_K)/norm(H(:,i_K));
[temp1,index1]=max(real(h_k'*C1));
h_RVQ = C1(:,index1)*norm(H(:,i_K));
QuantH(:,i_K)=h_RVQ;
end
err_RVQ=temp1;
else  fprintf('error');
err_RVQ = 0;
% err_RVQ = norm(H-QuantH,'fro')^2/norm(H,'fro')^2;
end

