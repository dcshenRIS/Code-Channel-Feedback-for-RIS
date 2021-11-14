function [H_SC,err_SC]= quantiz_SC_real(H,A,B)
%% descriptions of this function
% This function is to acquire the quan. channel accoring to the proposed scheme.
% ---------------- input descriptions -------------------------------------
%   "H" is the perfect channel with M(BS antennas) and K(users). 
%   "A" is the steering matrix.
%   "B" is the codebook size
% ---------------- output descriptions ------------------------------------
%   "H_SC" is is the fed back channel with codebook.
if B<=30
M = size(H, 1);
K = size(H, 2);
P = size(A, 2);
H_SC = zeros(size(H));
% Codebook
for i_K=1:1:K
C20=randn(P,2^B) + 1j*randn(P,2^B);
C2 = A(:,:,i_K) * C20;
C2 = C2./repmat(sqrt(sum(abs(C2).^2,1)),M,1);  
h_k = H(:,i_K)/norm(H(:,i_K));
[temp1,index1]=max(real(h_k'*C2));
method1=h_k'*C2(:,index1);
h_SC = C2(:,index1)*norm(H(:,i_K));
H_SC(:,i_K)=h_SC;
end
err_SC = temp1;
else  fprintf('error');
err_SC = 1;
end

