function [H_fb_SC_B0] = For_Parfor_Calcu_ProposedScheme_with_B0bits_quan_AoA_Rate(L1,on_gird_Hss_rum_trans,A_quan,B,N,M_long,K,spt_col,UM_long)
%% Calcu. the per-user rate for proposed scheme with quan. cascaded AoA and on-grid AoD(with errors)      
        for i_L1=1:L1
            H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
            A_B0(:,:,:)=A_quan(:,:,i_L1,:);%在这里区分的A_quan和A_total
            [H_SC_B0(:,i_L1,:),err_SC]= quantiz_SC_real(H,A_B0,B);%这里的A是B0量化后的
            %norm(permute(H_SC_B0(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
        end
        H_fb_SC_B0_rum=zeros(N,M_long,K,'double');
        
        H_fb_SC_B0_rum(:,spt_col,:)=H_SC_B0(:,:,:);
        for i_K=1:K
            H_fb_SC_B0(:,:,i_K)=H_fb_SC_B0_rum(:,:,i_K)*UM_long';
            
            %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
        end
       

end

