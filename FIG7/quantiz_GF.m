function [QuantH,err_GF]= quantiz_GF(H,A_quan,B,G)
if B<=20
M = size(H, 1);
K = size(H, 2);
P = size(G,1);
QuantH=zeros(size(H));
for i_K=1:1:K
C3=randn(P,2^B) + 1j*randn(P,2^B);
C3 = C3./repmat(sqrt(sum(abs(C3).^2,1)),P,1);  
g_k = G(:,i_K)/norm(G(:,i_K));
[temp3,index3]=max(real(g_k'*C3));
g_k_quan = C3(:,index3)*norm(G(:,i_K));
QuantH(:,i_K)= A_quan(:,:,i_K)*g_k_quan;
end
err_GF=temp3;
else  fprintf('error');
err_GF = 0;
% err_RVQ = norm(H-QuantH,'fro')^2/norm(H,'fro')^2;
end
