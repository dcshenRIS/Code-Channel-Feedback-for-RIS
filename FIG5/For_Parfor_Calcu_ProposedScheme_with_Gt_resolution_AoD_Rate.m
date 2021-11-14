function [H_fb_quan_on_gird] = For_Parfor_Calcu_ProposedScheme_with_Gt_resolution_AoD_Rate(L1,quan_on_gird_Hss_rum_trans,A_total,B,N,Value_M_long_table_of_i_M_long_index,K,spt_col_with_Mvar,UM_with_M_var)
       for i_L1=1:L1
            H(:,:)=quan_on_gird_Hss_rum_trans(:,i_L1,:);
            A(:,:,:)=A_total(:,:,i_L1,:);
            [H_SC_quan_on_gird(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
        end
        H_fb_SC_rum_quan_on_gird=zeros(N,Value_M_long_table_of_i_M_long_index,K);
        H_fb_SC_rum_quan_on_gird(:,spt_col_with_Mvar,:)=H_SC_quan_on_gird(:,:,:);
        

        
        for i_K=1:K
            H_fb_quan_on_gird(:,:,i_K)=H_fb_SC_rum_quan_on_gird(:,:,i_K)*UM_with_M_var';
        end
       

       

end

