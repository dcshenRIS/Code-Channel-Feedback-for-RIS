clear all;
%% Para. for Monte Carlo simulation
loop_num =250;

%% System Para.
M=64;% Num of antennas at BS
N=256;% Num of RIS elements
N1=16;
N2=N/N1;
L1=4;% Num of BS-RIS paths
L2=2;% Num of RIS-UE paths
Lr=4;% Num of BS-UE paths
K=4;% Num of UEs
P=L2;
d=0.5;

%% Para. setting for channel feedback
B0=6;%
%B0=12;%Quan. bits for cascaded AoAs
M_long=1024;%AoD resolution at the BS
%B_table = 5:2:15;%Codebook size
% B_table = 1:2:11;%Codebook size

%% Para. for CEO
N_iter=1000;

% %% Init.
% Rs0_table=zeros(loop_num,length(B_table));
% Rs_SC_table=zeros(loop_num,length(B_table));
% Rs_CRVQ_table=zeros(loop_num,length(B_table));
% Rs_CRVQ1_table=zeros(loop_num,length(B_table));

%% Gene dictionary matrix
UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFTæÿ’Û M*M_long
UM_long=(1/sqrt(M))*UM_long;

% Length_of_B_table=length(B_table);


SNR_table = -10:3:12;%SNR
Length_of_SNR_table=length(SNR_table);
%% Init.
Rs0_table=zeros(loop_num,length(SNR_table));
Rs_SC_table_with_1stCodebookSize=zeros(loop_num,length(SNR_table));
Rs_SC_table_with_2ndCodebookSize=zeros(loop_num,length(SNR_table));
Rs_SC_table_with_nofeedback=zeros(loop_num,length(SNR_table));



parfor i_loop=1:loop_num
    %% Gene. cascaded channel and direct channel
    [H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%H_d M*K
    [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_long);%B0_table(i_B0),M_long);
    %% Perfect total channel
    [~,~,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);
    H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA);

    
    for i_SNR=1:1:Length_of_SNR_table
        %% SNR
        norm_sq_sum=0;
        for i_K=1:1:K
            norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
        end
        
        norm_h_sq=norm_sq_sum/K;
        SNR = SNR_table(i_SNR);
        snr=10^(SNR/10);
        Gamma = snr*K/norm_h_sq;
        %B = SNR+5;
        p=Gamma/K;
        
        %% Calcu. the per-user rate with perfect CSIT
        Rs0_table(i_loop,i_SNR)=Rs_cal(H_pf_total,H_pf_total,p);
        
        %% Calcu. the per-user rate for proposed scheme with codebooksize B=4  (Total T/2 bits overhead in Fig 7)
%         Value_B_table_i_B=B_table(i_SNR);
        Value_B_table_1st_B=4;
        [H_fb_SC] = For_Parfor_calcu_proposed_scheme(L1,on_gird_Hss_rum_trans,A_quan,Value_B_table_1st_B,N,M_long,K,spt_col,UM_long)
        [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table_with_1stCodebookSize(i_loop,i_SNR)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
        
        
       %% Calcu. the per-user rate for proposed scheme with codebooksize B=8 (Total T bits overhead in Fig 7)
%         Value_B_table_i_B=B_table(i_SNR);
        Value_B_table_2nd_B=8;
        [H_fb_SC] = For_Parfor_calcu_proposed_scheme(L1,on_gird_Hss_rum_trans,A_quan,Value_B_table_2nd_B,N,M_long,K,spt_col,UM_long)
        [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table_with_2ndCodebookSize(i_loop,i_SNR)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
     
        
       %% Calcu. the per-user rate for proposed scheme with no feedback
%         Value_B_table_i_B=B_table(i_SNR);
       % Value_B_table_2nd_B=8;
        [H_fb_SC] = For_Parfor_calcu_proposed_scheme(L1,on_gird_Hss_rum_trans,A_quan,Value_B_table_2nd_B,N,M_long,K,spt_col,UM_long)
        rfactor=randn(size(H_fb_SC));%Random on behalf of  no feedback 
        [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(rfactor.*H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table_with_nofeedback(i_loop,i_SNR)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
        
  
      


    
 
   
     
        


        %% Show state.
        fprintf('For FIG7:i_loop=%d of %d,i_SNR =%d\n',i_loop,loop_num,SNR_table(i_SNR));
    end
end
%% Calcu. the average result with Monte Carlo
mRs0=mean(Rs0_table,1);
mRs_SC_with_1stCodebookSize=mean(Rs_SC_table_with_1stCodebookSize,1);
mRs_SC_with_2ndCodebookSize=mean(Rs_SC_table_with_2ndCodebookSize,1);
mRs_SC_with_nofeedback=mean(Rs_SC_table_with_nofeedback,1);
%% Plot
figure;
plot(SNR_table,mRs0,'k--','LineWidth',1.5);%14 bits means that we choose AoD resolution as 1024 and AoA quan. as 8 bits, then we can calcu. the 14 according to the formula for this paper.
hold on;
plot(SNR_table,mRs_SC_with_1stCodebookSize,'bd-','LineWidth',1.5);
hold on;
plot(SNR_table,mRs_SC_with_2ndCodebookSize,'rd-','LineWidth',1.5);
hold on;


plot(SNR_table,mRs_SC_with_nofeedback,'g--','LineWidth',1.5);
hold on;
legend('Perfect CSIT','Proposed scheme with T bits overhead','Proposed scheme with T/2 bits overhead','Without channel feedback');
% axis([3 15 0 4]);
grid on;
xlabel('SNR');
ylabel('Achievable sum-rate [bps/Hz]');



