clear all;
%% Para. for Monte Carlo simulation times
loop_num =1500;

%% System Para.
M=64;% Num of antennas at BS ULA
N=256;% Num of RIS elements UPA
N1=16;
N2=N/N1;
L1=4;% Num of BS-RIS paths
L2=2;% Num of RIS-UE paths
Lr=4;% Num of BS-UE paths
K=4;% Num of UEs
P=L2;
d=0.5;%distance of antenna 

%% Para. setting for channel feedback
B0=8;%Quan. bits for cascaded AoAs 
%B0=12;%Quan. bits for cascaded AoAs
M_long=1024;%AoD resolution at the BS
%B_table = 5:2:15;%Codebook size
B_table = 1:2:11;%Codebook size

%% Para. for CEO Algorithm
N_iter=1000;

%% Init.
Rs0_table=zeros(loop_num,length(B_table));
Rs_SC_table=zeros(loop_num,length(B_table));
Rs_CRVQ_table=zeros(loop_num,length(B_table));
Rs_CRVQ1_table=zeros(loop_num,length(B_table));

%% Gene dictionary matrix
UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFT Matrix with M_long resolution M*M_long
UM_long=(1/sqrt(M))*UM_long;

Length_of_B_table=length(B_table);
parfor i_loop=1:loop_num
    for i_B=1:1:Length_of_B_table
        %% Gene. cascaded channel and direct channel
        [H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%Direct channels H_d M*K
        [Hss,on_gird_Hss_rum_trans,spt_col,A_total,A_quan,g_total]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_long);%HSS cascaded channels N*M*K 
        %% Perfect total channel
        [~,~,THETA]=ACE_RIS_2bit(N_iter,Hss,H_d);%Passive beamforming for RIS
        H_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA); %Total channel (perfectCSI)
        norm_sq_sum=0;
        for i_K=1:1:K
            norm_sq_sum = norm_sq_sum + norm(H_pf_total(:,i_K))^2;
        end
        
        %% Calculate for transmit energy 
        norm_h_sq=norm_sq_sum/K;%Channel Ener.
        SNR = 5;
        snr=10^(SNR/10);%Receive SNR (NoisePower=1)
        Gamma = snr*K/norm_h_sq;%TransmitEnergy
        %B = SNR+5;
        p=Gamma/K;%TransmitEnergy for per user
        
        %% Calcu. the per-user rate with perfect CSIT
        Rs0_table(i_loop,i_B)=Rs_cal(H_pf_total,H_pf_total,p);%
        %% Calcu. the per-user rate for proposed scheme
        Value_B_table_i_B=B_table(i_B);
        [H_fb_SC] = For_Parfor_calcu_proposed_scheme(L1,on_gird_Hss_rum_trans,A_quan,Value_B_table_i_B,N,M_long,K,spt_col,UM_long)
        [~,D,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);%Passive BF according to the feedback SCI
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table(i_loop,i_B)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
        
        %% Calcu. the per-user rate for Statistics-based codebook feedback(The blue curve in FIG4)
        H_rum=zeros(N,M_long,K);
        H_rum(:,spt_col,:)=on_gird_Hss_rum_trans(:,:,:);
        [H_CRVQ,H_CRVQ1] = For_Parfor_calcu_baseline_scheme(K,H_rum,UM_long,N);%The Statistics-based scheme with two different codebook size 
        
        [~,D,THETA_CRVQ]=ACE_RIS_2bit(N_iter,H_CRVQ,H_d);
        [~,D1,THETA_CRVQ1]=ACE_RIS_2bit(N_iter,H_CRVQ1,H_d);
        H_fb_CRVQ_total = RIS_Change_Cascaded_to_Total_channel(H_CRVQ,H_d,THETA_CRVQ);
        H_fb_CRVQ_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_CRVQ);
        
        H_fb_CRVQ1_total = RIS_Change_Cascaded_to_Total_channel(H_CRVQ1,H_d,THETA_CRVQ1);
        H_fb_CRVQ1_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_CRVQ1);
        
        Rs_CRVQ_table(i_loop,i_B)=Rs_cal(H_fb_CRVQ_pf_total,H_fb_CRVQ_total,p);
        Rs_CRVQ1_table(i_loop,i_B)=Rs_cal(H_fb_CRVQ1_pf_total,H_fb_CRVQ1_total,p);
        %% Show state.
        fprintf('For FIG4:i_loop=%d of %d,i_proposed_codebook_size B=%d\n',i_loop,loop_num,B_table(i_B));
        
    end
end
%% Save the simulation results
filename   =   strcat('Rs0_table',num2str(datestr(now,30)),   '.mat');
save(['.\Data\FIG4data\',filename],    'Rs0_table')%,'-v7.3');
filename   =   strcat('Rs_SC_table',num2str(datestr(now,30)),   '.mat');
save(['.\Data\FIG4data\',filename],    'Rs_SC_table')%,'-v7.3');
filename   =   strcat('Rs_CRVQ_table',num2str(datestr(now,30)),   '.mat');
save(['.\Data\FIG4data\',filename],    'Rs_CRVQ_table')%,'-v7.3');
filename   =   strcat('Rs_CRVQ1_table',num2str(datestr(now,30)),   '.mat');
save(['.\Data\FIG4data\',filename],    'Rs_CRVQ1_table')%,'-v7.3');

%% Calcu. the average result with Monte Carlo
mRs0=mean(Rs0_table,1);
mRs_SC=mean(Rs_SC_table,1);
mRs_CRVQ=mean(Rs_CRVQ_table,1);
mRs_CRVQ1=mean(Rs_CRVQ1_table,1);%320bits
%% Plot
figure;
plot(4*B_table+14,mRs0,'k--','LineWidth',1.5);%14 bits means that we choose AoD resolution as 1024 and AoA quan. as 8 bits, then we can calcu. the 14 according to the formula (25) for this paper.
hold on;
plot(4*B_table+14,mRs_SC,'rd-','LineWidth',1.5);
hold on;
plot(4*B_table+14,mRs_CRVQ1,'bo-','LineWidth',1.5);
hold on;
plot(4*B_table+14,mRs_CRVQ,'b^-','LineWidth',1.5);
hold on;


legend('Perfect CSIT','Proposed dimension reduced channel feedback scheme','Statistics-based channel feedback with overhead=320 bits','Statistics-based channel feedback with overhead=64 bits');
% axis([3 15 0 4]);
grid on;
xlabel('Per-user feedback overhead for proposed dimension reduced channel feedback (bits)');
ylabel('Achievable sum-rate [bps/Hz]');


