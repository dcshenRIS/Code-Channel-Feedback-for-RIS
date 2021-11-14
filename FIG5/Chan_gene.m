function [H,A,A_quan,G]=Chan_gene(M,K,P,B0)
%% descriptions of this function
% This function is for the direct channel generation.
% ---------------- input descriptions -------------------------------------
%   "M" is the number of antennas of BS.
%   "K" is the number of UEs.
%   "P" is the number of paths.
%   "B0" is the quan. bits. (but not used in our work)
% 
% ---------------- output descriptions ------------------------------------
%   "H" is the direct channel.
%% Gene channel
d = 0.5;
sin_theta_table = zeros(K,P);
sin_theta_table_quan = zeros(K,P);

for i_K=1:1:K
    theta = random('unif',-1/2,1/2,1,P)*pi;%continuous channel angles w/o girds

    sin_theta_table(i_K,:)=sin(theta);
    sin_theta_table_quan(i_K,:)=UQ(sin_theta_table(i_K,:),B0,-1,1);

end
H=zeros(M,K);
A = zeros(M,P,K);
A_quan = zeros(M,P,K);
G = zeros(P,K);
for i_K=1:1:K
        A(:,:,i_K) = sqrt(1/M)*exp([0:1:M-1]'*(1j*2*pi*d)*sin_theta_table(i_K,:));
        A_quan(:,:,i_K) = sqrt(1/M)*exp([0:1:M-1]'*(1j*2*pi*d)*sin_theta_table_quan(i_K,:));

        g = randn(P,1) + 1j*randn(P,1);
        G(:,i_K)=g;
        h  = A(:,:,i_K)*g;
        H(:,i_K) = h;  
end