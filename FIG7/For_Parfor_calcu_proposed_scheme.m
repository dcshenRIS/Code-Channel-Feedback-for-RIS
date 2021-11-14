function [H_fb_SC] = For_Parfor_calcu_proposed_scheme(L1,on_gird_Hss_rum_trans,A_quan,Value_B_table_i_B,N,M_long,K,spt_col,UM_long)
        %% Calcu. the per-user rate for proposed scheme 
        for i_L1=1:L1
            H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
            A(:,:,:)=A_quan(:,:,i_L1,:);
            [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,Value_B_table_i_B);
            %norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
        end
        H_fb_SC_rum=zeros(N,M_long,K);
        H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
        for i_K=1:K
            H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
            
            %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
        end
end

