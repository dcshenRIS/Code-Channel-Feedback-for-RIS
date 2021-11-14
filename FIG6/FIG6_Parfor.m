clear all;
%% Para. for Monte Carlo simulation
loop_num = 1000;

%% System Para.
M=64;
N=256;
N1=16;
N2=N/N1;
L1=4;
L2=2;
Lr=4;
K=4;
P=L2;
d=0.5;

%% Para. setting for channel feedback
M_long=1024;%AoD resolution at the BS
B=14;%Codebook size
B0_table = 1:1:8;   %Quan. bits for cascaded AoAs

%% Para. for CEO
N_iter=2000;



%% Init.
Rs0_table=zeros(loop_num,length(B0_table));
Rs_SC_table=zeros(loop_num,length(B0_table));
Rs_SC_B0_table=zeros(loop_num,length(B0_table));
H_SC=zeros(N,L1,K,'double');
H_SC_B0=zeros(N,L1,K,'double');
H_fb_SC=zeros(N,M,K,'double');
H_fb_SC_B0=zeros(N,M,K,'double');
%% Gene dictionary matrix
UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFTæÿ’Û M*M
UM_long=(1/sqrt(M))*UM_long;

Length_of_B0_table=length(B0_table);
parfor i_loop=1:loop_num
    for i_B0=1:1:Length_of_B0_table%length(B0_table)
        %% Gene. cascaded channel and direct channel
        [H_d,~,~,~]=Chan_gene(M,K,Lr,7);%B0_table(i_B0));%H_d M*K
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
        
        p=Gamma/K;
        
        
        
        %% Perfect Rate
        Rs0_table(i_loop,i_B0)=Rs_cal(H_pf_total,H_pf_total,p);
        %% Calcu. the per-user rate for proposed scheme with perfect cascaded AoA but on-grid AoD(with errors)
        [H_fb_SC] = For_Parfor_Calcu_ProposedScheme_with_perfect_AoA_Rate(L1,on_gird_Hss_rum_trans,A_total,B,N,M_long,K,spt_col,UM_long);
        
        [~,~,THETA_SC]=ACE_RIS_2bit(N_iter,H_fb_SC,H_d);
        H_fb_SC_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC,H_d,THETA_SC);
        H_fb_SC_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC);
        %Cal Rate
        Rs_SC_table(i_loop,i_B0)=Rs_cal(H_fb_SC_pf_total,H_fb_SC_total,p);
        %% Calcu. the per-user rate for proposed scheme with quan. cascaded AoA and on-grid AoD(with errors)
        [H_fb_SC_B0] = For_Parfor_Calcu_ProposedScheme_with_B0bits_quan_AoA_Rate(L1,on_gird_Hss_rum_trans,A_quan,B,N,M_long,K,spt_col,UM_long);
        
        [~,D,THETA_SC_B0]=ACE_RIS_2bit(N_iter,H_fb_SC_B0,H_d);
        H_fb_SC_B0_total = RIS_Change_Cascaded_to_Total_channel(H_fb_SC_B0,H_d,THETA_SC_B0);
        H_fb_SC_B0_pf_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_SC_B0);
        %Cal Rate
        Rs_SC_B0_table(i_loop,i_B0)=Rs_cal(H_fb_SC_B0_pf_total,H_fb_SC_B0_total,p);
        
        
        %% Show state.
        fprintf('For FIG6:i_loop=%d of %d,cascaded AoA quan. bit B0=%d\n',i_loop,loop_num,i_B0);
        
        
    end
end
%% Save the simulation results
filename   =   strcat('Rs0_table',   '.mat');
save(['.\Data\FIG6data\',filename],    'Rs0_table')%,'-v7.3');

filename   =   strcat('Rs_SC_table',   '.mat');
save(['.\Data\FIG6data\',filename],    'Rs_SC_table')%,'-v7.3');

filename   =   strcat('Rs_SC_B0_table',   '.mat');
save(['.\Data\FIG6data\',filename],    'Rs_SC_B0_table')%,'-v7.3');
%% Calcu. the average result with Monte Carlo
mRs0=mean(Rs0_table,1);
mRs_SC=mean(Rs_SC_table,1);
mRs_SC_B0=mean(Rs_SC_B0_table,1);
%  mRs0=mean(Rs0_table,1);
%  mRs_SC=mean(Rs_SC_table,1);
% mRs_SC_B0=mean(Rs_SC_B0_table,1);
%
%% Plot
figure;
plot(B0_table,mRs0,'k--','LineWidth',1.5);
hold on;
% hold on;
plot(B0_table,mRs_SC,'rd-','LineWidth',1.5);
hold on;
plot(B0_table,mRs_SC_B0,'b^-','LineWidth',1.5);
hold on;


legend('Perfect CSIT','Proposed scheme with  perfect cascaded AoAs ','Proposed scheme with  quantized cascaded AoAs');
% axis([3 15 0 4]);
grid on;
xlabel('Cascaded AoA quanzition bits B_0 at the RIS');
ylabel('Achievable sum-rate [bps/Hz]');

