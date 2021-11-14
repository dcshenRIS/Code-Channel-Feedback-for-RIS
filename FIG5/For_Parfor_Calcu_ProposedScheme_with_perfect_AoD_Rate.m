function [H_fb_quan] = For_Parfor_Calcu_ProposedScheme_with_perfect_AoD_Rate(L1,quan_Hss_rum_trans,A_total,B,N,M_long_max,K,spt_col_with_Mlong,UM_with_M_long)
for i_L1=1:L1
    H(:,:)=quan_Hss_rum_trans(:,i_L1,:);
    A(:,:,:)=A_total(:,:,i_L1,:);
    [H_SC_quan(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
end
H_fb_SC_rum_quan=zeros(N,M_long_max,K);
H_fb_SC_rum_quan(:,spt_col_with_Mlong,:)=H_SC_quan(:,:,:);



for i_K=1:K
    H_fb_quan(:,:,i_K)=H_fb_SC_rum_quan(:,:,i_K)*UM_with_M_long';
    
end
end

