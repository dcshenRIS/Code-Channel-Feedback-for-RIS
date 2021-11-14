function [H_CRVQ,H_CRVQ1] = For_Parfor_calcu_baseline_scheme(K,H_rum,UM_long,N)
for i_K=1:K
    Hss_og(:,:,i_K)=H_rum(:,:,i_K)*UM_long';
    
end
for i_N=1:N
    
    H1(:,:)=Hss_og(i_N,:,:);
    %% Channel Covariance
    R_t=H1 * H1';
    [U,S,V] = svd(R_t);
    R_sqt = U*sqrt(S)*U';
     [H_CRVQ(i_N,:,:),err_CRVQ]= quantiz_CRVQ_real(H1,R_sqt,1);% with codebook size 1bit per vector
     [H_CRVQ1(i_N,:,:),err_CRVQ1]= quantiz_CRVQ_real(H1,R_sqt,5);% with codebook size 5bit per vector
end
end

