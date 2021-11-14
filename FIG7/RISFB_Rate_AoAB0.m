clear all;
M=64;
N=256;
N1=16;
N2=N/N1;
L1=4;
L2=2;
Lr=4;
K=4;
%B0=7;
M_long=1024;%512;
P=L2;
B=10;%实际上没有用
d=0.5;
N_iter=1000;

loop_num = 30;
B0_table = 1:1:8;   

Rs0_table=zeros(loop_num,length(B0_table));
Rs_SC_table=zeros(loop_num,length(B0_table));
Rs_SC_B0_table=zeros(loop_num,length(B0_table));
%% Gene DFT Channel (Right DFT UM)

UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFT矩阵 M*M
UM_long=(1/sqrt(M))*UM_long;
[H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%H_d M*K

H_SC=zeros(N,L1,K,'double');
H_SC_B0=zeros(N,L1,K,'double');
H_fb_SC=zeros(N,M,K,'double');
H_fb_SC_B0=zeros(N,M,K,'double');

for i_loop=1:loop_num
    for i_B0=1:1:length(B0_table)
        [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0_table(i_B0),M_long);%B0_table(i_B0),M_long);
        %% Perfect total channel
        [~,~,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
        H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);
        norm_sq_sum=0;
        for i_K=1:1:K
            norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
        end

        norm_h_sq=norm_sq_sum/K;
        SNR = 5;
        snr=10^(SNR/10);
        Gamma = snr*K/norm_h_sq;
        B = SNR+6;
        p=Gamma/K;
                %% Perfect Rate
        Rs0_table(i_loop,i_B0)=Rs_cal(H_pf_total,H_pf_total,p);
        %Rs0_table(i_loop,i_B0)=Rs0_table(i_loop,1);
        %% SC with Perfect AoD Rate
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
        [~,~,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
        NMSE1=10*log10(norm(H_fb_SC_pf_total-H_fb_SC_total,'fro')^2/norm(H_fb_SC_pf_total,'fro')^2);
        %Rs_SC_table(i_loop,i_B0)=Rs_SC_table(i_loop,1);
         %% SC with Quan B0 AoD Rate       
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
        [~,D,THETA_SC_B0]=ACE_RIS_2bit(N_iter,H_fb_SC_B0,H_d);
        H_fb_SC_B0_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC_B0,H_d,THETA_SC_B0);
        H_fb_SC_B0_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC_B0);
        %Cal Rate
        Rs_SC_B0_table(i_loop,i_B0)=Rs_cal(H_fb_SC_B0_pf_total,H_fb_SC_B0_total,p);
        
     fprintf('i_loop=%d,AoA quan B0=%d\n',i_loop,i_B0); 

    
    end
%     for i_B0=1:1:length(B0_table)
%         Rs0_table(i_loop,i_B0)=Rs0_table(i_loop,length(B0_table));
%         Rs_SC_table(i_loop,i_B0)=Rs_SC_table(i_loop,length(B0_table));
%     end
end
 mRs0=mean(Rs0_table,1)/K;
 mRs_SC=mean(Rs_SC_table,1)/K;
mRs_SC_B0=mean(Rs_SC_B0_table,1)/K;
% 
figure;
plot(B0_table,mRs0,'kp--','LineWidth',1.5);
hold on;
% hold on;
plot(B0_table,mRs_SC,'rd-','LineWidth',1.5);
hold on;
plot(B0_table,mRs_SC_B0,'c^-','LineWidth',1.5);
hold on;


legend('Perfect CSIT','Proposed scheme with  perfect cascaded AoAs ','Proposed scheme with  imperfect cascaded AoAs');
% axis([3 15 0 4]);
grid on;
xlabel('Cascaded AoA quanzition bits B_0 at the RIS');
ylabel('Per-user rate (bps/Hz)');


% for i_loop=1:loop_num
%     
%         [H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%H_d M*K
%         [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene_forRate_AoDB0(M,N,N1,N2,L1,L2,K,7,M_long);%B0_table(i_B0),M_long);
% for i_B0=1:1:length(B0_table)
% %                  [H_d,~,~,~]=Chan_gene(M,K,Lr,B0_table(i_B0));%H_d M*K
% %          [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene_forRate_AoDB0(M,N,N1,N2,L1,L2,K,B0_table(i_B0),M_long);
%         %% Perfect total channel
%         [~,D,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
% %         N_gird=N;
% %         theta=zeros(1,N);
% %         theta_gird=[0:1:(N_gird-1)]*(1/N_gird)*pi/2;
% %         theta_rand=randi(N_gird,[1,N]);
% %         THETA=exp(1j*theta_gird(theta_rand));%1*N
% %         %THETA=exp(1j*zeros(1,N));%1*N
% 
%         H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);
%         norm_sq_sum=0;
%         for i_K=1:1:K
%             norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
%         end
%         norm_h_sq=norm_sq_sum/K;
%         
%         
%         SNR = 10;
%         snr=10^(SNR/10);
%         Gamma = snr*K/norm_h_sq;
%         B = SNR;
%         p=Gamma/K;
%         
%         
% 
%         %% Perfect Rate
%         Rs0_table(i_loop,i_B0)=Rs_cal(H_pf_total,H_pf_total,p);
% %         %% SC with Perfect AoD Rate       
% %         for i_L1=1:L1
% %             H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %             A(:,:,:)=A_total(:,:,i_L1,:);
% %             [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
% %             norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2
% %         end
% %         H_fb_SC_rum=zeros(N,M_long,K);
% %         H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
% %         for i_K=1:K
% %             H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
% %             
% %             %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %         end
% %         [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
% %         H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
% %         H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
% %         %Cal Rate
% %         Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
% %         %% SC with Quan B0 AoD Rate       
% %         for i_L1=1:L1
% %             H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %             A_B0(:,:,:)=A_quan(:,:,i_L1,:);%在这里区分的A_quan和A_total
% %             [H_SC_B0(:,i_L1,:),err_SC]= quantiz_SC_real(H,A_B0,B);%这里的A是B0量化后的
% %             norm(permute(H_SC_B0(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% %         end
% %         H_fb_SC_B0_rum=zeros(N,M_long,K);
% %         H_fb_SC_B0_rum(:,spt_col,:)=H_SC_B0(:,:,:);
% %         for i_K=1:K
% %             H_fb_SC_B0(:,:,i_K)=H_fb_SC_B0_rum(:,:,i_K)*UM_long';
% %             
% %             %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %         end
% %         [~,D,THETA_SC_B0]=ACE_RIS_2bit(N_iter,H_fb_SC_B0,H_d);
% %         H_fb_SC_B0_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC_B0,H_d,THETA_SC_B0);
% %         H_fb_SC_B0_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC_B0);
% %         %Cal Rate
% %         Rs_SC_B0_table(i_loop,i_B0)=Rs_cal(H_fb_SC_B0_pf_total,H_fb_SC_B0_total,p);
%         
%      fprintf('i_loop=%d,B0=%d\n',i_loop,i_B0);   
% end
% end
% 
% 
% 
% 
% 
% % 
% %         N_gird=N;
% %         theta=zeros(1,N);
% %         theta_gird=[0:1:(N_gird-1)]*(1/N_gird)*pi/2;
% %         theta_rand=randi(N_gird,[1,N]);
% %         THETA=exp(1j*theta_gird(theta_rand));%1*N
% 
% % for i_B0=1:1:length(B0_table)
% %     [H_d,~,~,~]=Chan_gene(M,K,Lr,B0_table(i_B0));%B0_table(i_B0));%H_d M*K
% %     [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0_table(i_B0),M_long);%B0_table(i_B0),M_long);
% %     [~,D,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
% %     
% %         %% Perfect total channel        
% %          H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);
% %         norm_sq_sum=0;
% %         for i_K=1:1:K
% %             norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
% %         end
% %         norm_h_sq=norm_sq_sum/K;        
% %         SNR = 10;
% %         snr=10^(SNR/10);
% %         Gamma = snr*K/norm_h_sq;
% %         B = SNR;
% %         p=Gamma/K;
% %         for i_loop=1:loop_num
% %             %% Perfect Rate
% %             Rs0_table(i_loop,i_B0)=Rs_cal(H_pf_total,H_pf_total,p);
% %             %% SC with Perfect AoD Rate
% %             for i_L1=1:L1
% %                 H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %                 A(:,:,:)=A_total(:,:,i_L1,:);
% %                 [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
% %                 norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% %             end
% %             H_fb_SC_rum=zeros(N,M_long,K);
% %             H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
% %             for i_K=1:K
% %                 H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
% %                 
% %                 %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %             end
% %             [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
% %             H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
% %             H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
% %             %Cal Rate
% %             Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
% %             %% SC with Quan B0 AoD Rate
% %             for i_L1=1:L1
% %                 H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %                 A_B0(:,:,:)=A_quan(:,:,i_L1,:);%在这里区分的A_quan和A_total
% %                 [H_SC_B0(:,i_L1,:),err_SC]= quantiz_SC_real(H,A_B0,B);%这里的A是B0量化后的
% %                 norm(permute(H_SC_B0(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% %             end
% %             H_fb_SC_B0_rum=zeros(N,M_long,K);
% %             H_fb_SC_B0_rum(:,spt_col,:)=H_SC_B0(:,:,:);
% %             for i_K=1:K
% %                 H_fb_SC_B0(:,:,i_K)=H_fb_SC_B0_rum(:,:,i_K)*UM_long';
% %                 
% %                 %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %             end
% %             [~,D,THETA_SC_B0]=ACE_RIS_2bit(N_iter,H_fb_SC_B0,H_d);
% %             H_fb_SC_B0_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC_B0,H_d,THETA_SC_B0);
% %             H_fb_SC_B0_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC_B0);
% %             %Cal Rate
% %             Rs_SC_B0_table(i_loop,i_B0)=Rs_cal(H_fb_SC_B0_pf_total,H_fb_SC_B0_total,p);
% %             fprintf('i_loop=%d,B0=%d\n',i_loop,i_B0);
% %         end
% % 
% % 
% % end
% mRs0=mean(Rs0_table,1)/K;
% mRs_SC=mean(Rs_SC_table,1)/K;
% mRs_SC_B0=mean(Rs_SC_B0_table,1)/K;
% % 
% figure;
% plot(B0_table,mRs0,'kp--','LineWidth',1.5);
% hold on;
% hold on;
% plot(B0_table,mRs_SC,'rd-','LineWidth',1.5);
% hold on;
% plot(B0_table,mRs_SC_B0,'c^-','LineWidth',1.5);
% 
% 
% legend('Perfect CSIT','Proposed feedback strategy ','Conventional channel statistics-based codebook feedback','RVQ codebook feedback');
% % axis([3 15 0 4]);
% grid on;
% xlabel('SNR (dB)');
% ylabel('Per-user rate (bps/Hz)');
% 
% 
% 
% % 
% % 
% % 
% % for i_loop=1:loop_num
% %     
% %         [H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%H_d M*K
% %         [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene_forRate_AoDB0(M,N,N1,N2,L1,L2,K,7,M_long);%B0_table(i_B0),M_long);
% % for i_B0=1:1:length(B0_table)
% % %                  [H_d,~,~,~]=Chan_gene(M,K,Lr,B0_table(i_B0));%H_d M*K
% % %          [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene_forRate_AoDB0(M,N,N1,N2,L1,L2,K,B0_table(i_B0),M_long);
% %         %% Perfect total channel
% %         [~,D,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
% % %         N_gird=N;
% % %         theta=zeros(1,N);
% % %         theta_gird=[0:1:(N_gird-1)]*(1/N_gird)*pi/2;
% % %         theta_rand=randi(N_gird,[1,N]);
% % %         THETA=exp(1j*theta_gird(theta_rand));%1*N
% % %         %THETA=exp(1j*zeros(1,N));%1*N
% % 
% %         H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);
% %         norm_sq_sum=0;
% %         for i_K=1:1:K
% %             norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
% %         end
% %         norm_h_sq=norm_sq_sum/K;
% %         
% %         
% %         SNR = 10;
% %         snr=10^(SNR/10);
% %         Gamma = snr*K/norm_h_sq;
% %         B = SNR;
% %         p=Gamma/K;
% %         
% %         
% % 
% %         %% Perfect Rate
% %         Rs0_table(i_loop,i_B0)=Rs_cal(H_pf_total,H_pf_total,p);
% % %         %% SC with Perfect AoD Rate       
% % %         for i_L1=1:L1
% % %             H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% % %             A(:,:,:)=A_total(:,:,i_L1,:);
% % %             [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
% % %             norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2
% % %         end
% % %         H_fb_SC_rum=zeros(N,M_long,K);
% % %         H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
% % %         for i_K=1:K
% % %             H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
% % %             
% % %             %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% % %         end
% % %         [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
% % %         H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
% % %         H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
% % %         %Cal Rate
% % %         Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
% % %         %% SC with Quan B0 AoD Rate       
% % %         for i_L1=1:L1
% % %             H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% % %             A_B0(:,:,:)=A_quan(:,:,i_L1,:);%在这里区分的A_quan和A_total
% % %             [H_SC_B0(:,i_L1,:),err_SC]= quantiz_SC_real(H,A_B0,B);%这里的A是B0量化后的
% % %             norm(permute(H_SC_B0(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% % %         end
% % %         H_fb_SC_B0_rum=zeros(N,M_long,K);
% % %         H_fb_SC_B0_rum(:,spt_col,:)=H_SC_B0(:,:,:);
% % %         for i_K=1:K
% % %             H_fb_SC_B0(:,:,i_K)=H_fb_SC_B0_rum(:,:,i_K)*UM_long';
% % %             
% % %             %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% % %         end
% % %         [~,D,THETA_SC_B0]=ACE_RIS_2bit(N_iter,H_fb_SC_B0,H_d);
% % %         H_fb_SC_B0_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC_B0,H_d,THETA_SC_B0);
% % %         H_fb_SC_B0_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC_B0);
% % %         %Cal Rate
% % %         Rs_SC_B0_table(i_loop,i_B0)=Rs_cal(H_fb_SC_B0_pf_total,H_fb_SC_B0_total,p);
% %         
% %      fprintf('i_loop=%d,B0=%d\n',i_loop,i_B0);   
% % end
% % end
% 
% 
% 
% % SNR_table = 1:3:16;
% % Rs0_table=zeros(loop_num,length(SNR_table));
% % Rs_SC_table=zeros(loop_num,length(SNR_table));
% % Rs_RVQ_table=zeros(loop_num,length(SNR_table));
% % 
% % for i_loop=1:loop_num
% %     [H_d,~,~,~]=Chan_gene(M,K,Lr,B0);%H_d M*K
% %     [Hss,on_gird_Hss_rum_trans,spt_col,A_total,g_total]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_long);
% %     %% Gene DFT Channel (Right DFT UM)
% %     
% %     UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFT矩阵 M*M
% %     UM_long=(1/sqrt(M))*UM_long;
% %     
% %     %% Perfect total channel
% %     [~,D,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
% %     H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);
% %     norm_sq_sum=0;
% %     for i_K=1:1:K
% %         norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
% %     end
% %     norm_h_sq=norm_sq_sum/K;
% %     
% %     
% %     for i_SNR=1:1:length(SNR_table)
% %     SNR = SNR_table(i_SNR);
% %     snr=10^(SNR/10);
% %     Gamma = snr*K/norm_h_sq;
% %     B = SNR+5;
% %     B_table(i_SNR) = B;
% %     %% SC method
% %     for i_L1=1:L1
% %         H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %         A(:,:,:)=A_total(:,:,i_L1,:);
% %         [H_SC(:,i_L1,:),err_SC]= quantiz_SC_real(H,A,B);
% %         norm(permute(H_SC(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2
% %     end
% %     H_fb_SC_rum=zeros(N,M_long,K);
% %     H_fb_SC_rum(:,spt_col,:)=H_SC(:,:,:);
% %     for i_K=1:K
% %         H_fb_SC(:,:,i_K)=H_fb_SC_rum(:,:,i_K)*UM_long';
% %         
% %         %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %     end
% %     [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
% %     H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
% %     H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
% %     %% RVQ method
% %     for i_L1=1:L1
% %         H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %         [H_RVQ(:,i_L1,:),err_RVQ]= quantiz_RVQ_real(H, B);
% %         norm(permute(H_RVQ(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% %     end
% %     H_fb_RVQ_rum=zeros(N,M_long,K);
% %     H_fb_RVQ_rum(:,spt_col,:)=H_RVQ(:,:,:);
% %     for i_K=1:K
% %         H_fb_RVQ(:,:,i_K)=H_fb_RVQ_rum(:,:,i_K)*UM_long';
% %         
% %         %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %     end
% %     [~,D,THETA_RVQ]=ACE_RIS_2bit(N_iter,H_fb_RVQ,H_d);
% %     H_fb_RVQ_total = RIS_Change_Cascaded_to_Total_channel(H_fb_RVQ,H_d,THETA_RVQ);
% %     H_fb_RVQ_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_RVQ);
% %     
% %     %% CRVQ method
% %     for i_L1=1:L1
% %         H(:,:)=on_gird_Hss_rum_trans(:,i_L1,:);
% %         %% Channel Covariance
% %         R_t=H * H';
% %         [U,S,V] = svd(R_t);
% %         R_sqt = U*sqrt(S)*U';
% %         
% %         [H_CRVQ(:,i_L1,:),err_CRVQ]= quantiz_CRVQ_real(H,R_sqt,B);
% %         norm(permute(H_CRVQ(:,i_L1,:),[1,3,2])-H,'fro')^2/norm(H,'fro')^2;
% %     end
% %     H_fb_CRVQ_rum=zeros(N,M_long,K);
% %     H_fb_CRVQ_rum(:,spt_col,:)=H_CRVQ(:,:,:);
% %     for i_K=1:K
% %         H_fb_CRVQ(:,:,i_K)=H_fb_CRVQ_rum(:,:,i_K)*UM_long';
% %         
% %         %     norm(H_fb_SC(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2
% %     end
% %     [~,D,THETA_CRVQ]=ACE_RIS_2bit(N_iter,H_fb_CRVQ,H_d);
% %     H_fb_CRVQ_total = RIS_Change_Cascaded_to_Total_channel(H_fb_CRVQ,H_d,THETA_CRVQ);
% %     H_fb_CRVQ_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_CRVQ);
% %     
% %      p=Gamma/K;
% % Rs0_table(i_loop,i_SNR)=Rs_cal(H_pf_total,H_pf_total,p);
% % Rs_SC_table(i_loop,i_SNR)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
% % Rs_RVQ_table(i_loop,i_SNR)=Rs_cal(H_fb_RVQ_pf_total,H_fb_RVQ_total,p);
% % Rs_CRVQ_table(i_loop,i_SNR)=Rs_cal(H_fb_CRVQ_pf_total,H_fb_CRVQ_total,p);
% %     end
% % end
% % mRs0=mean(Rs0_table,1)/K;
% % mRs_SC=mean(Rs_SC_table,1)/K;
% % mRs_CRVQ=mean(Rs_CRVQ_table,1)/K;
% % mRs_RVQ=mean(Rs_RVQ_table,1)/K;
% 
