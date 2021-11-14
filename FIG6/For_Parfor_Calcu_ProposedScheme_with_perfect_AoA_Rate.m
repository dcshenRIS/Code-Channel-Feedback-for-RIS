function [H_fb_SC] = For_Parfor_Calcu_ProposedScheme_with_perfect_AoA_Rate(L1,on_gird_Hss_rum_trans,A_total,B,N,M_long,K,spt_col,UM_long)
        for i_L1=1:L1
            H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
            A(:,:,:)=A_total(:,:,i_L1,:);
            [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
            %norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
        end
        H_fb_SC_rum=zeros(N,M_long,K,'double');
        
        H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
        for i_K=1:K
            H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
            
            %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
        end
%         [~,~,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
%         H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
%         H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
%         %Cal Rate
%         Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
end

