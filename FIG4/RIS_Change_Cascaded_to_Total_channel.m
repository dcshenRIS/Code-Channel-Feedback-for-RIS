function [H_total] = RIS_Change_Cascaded_to_Total_channel(H_refl,H_d,THETA)
%% descriptions of this function
% This function is to acquire the total channel between the BS and UE.
% ---------------- input descriptions -------------------------------------
%   "H_refl" is the cascaded channel. 
%   "H_d" is the direct channel.
%   "THETA" is the phase shift of RIS.
% ---------------- output descriptions ------------------------------------
%   "H_total" is is the total channel between the BS and UE.
N=size(H_refl,1);
M=size(H_refl,2);
K=size(H_refl,3);

for i_K=1:K
    H_pf_refl_link(:,i_K)=THETA*H_refl(:,:,i_K);
end
H_total=H_pf_refl_link+H_d;
end

