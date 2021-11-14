clear all;

%% Monte Carlo times
loop_num =1600;

%% Para. setting
M=64; %Num of TX ports
N=256; %Num of RIS elements
N1=16;
N2=N/N1;

L1=4;  %Path num of BS-RIS
L2=2;  %Path num of UE-RIS
Lr=4;  %Path num of BS-UE

K=4;  %Num of UE
d=0.5;
%% For CEO
N_iter=1000;
%% FB Para.
B0=8;%cascaded AoA quan. bits
%P=L2;
B=11; %Codebook size
M_long_max=2048;%Max resolution for AoD
M_long_index_table = log2(16):1:log2(M_long_max);
%M_index_table=M_long_index_table(1:length(M_long_index_table)-1);
M_long_table=2.^M_long_index_table;


%% Per-user rate
Rs0_table=zeros(loop_num,length(M_long_index_table));
Rs1_table=zeros(loop_num,length(M_long_index_table));
Rs2_table=zeros(loop_num,length(M_long_index_table));




UM_with_M_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long_max-1)/2:1:(M_long_max/2)]*(2/M_long_max));       %DFT矩阵 M*M
UM_with_M_long=(1/sqrt(M))*UM_with_M_long;


H_SC=zeros(N,L1,K,'double');
H_SC_M_long=zeros(N,L1,K,'double');
H_fb_SC=zeros(N,M,K,'double');
on_gird_Hss_pf_fb=zeros(N,M,K,'double');
H_fb_SC_M_long=zeros(N,M,K,'double');

Length_M_long_index_table=length(M_long_index_table);
parfor i_loop=1:loop_num
    for i_M_long_index=1:1:Length_M_long_index_table
        
        UM_with_M_var=exp(1i*2*pi*[0:M-1]'*d*[-(M_long_table(i_M_long_index)-1)/2:1:(M_long_table(i_M_long_index)/2)]*(2/M_long_table(i_M_long_index)));       %DFT矩阵 M*M_long
        UM_with_M_var=(1/sqrt(M))*UM_with_M_var;
        
        %% Gene direct Channel
        [H_d,~,~,~]=Chan_gene(M,K,Lr,B0); %Gnen the direct channel %H_d M*K
        [Hss,quan_Hss_rum_trans,quan_on_gird_Hss_rum_trans,A_total,spt_col_with_Mlong,spt_col_with_Mvar]=RIS_New_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_long_table(i_M_long_index),M_long_max);
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
        Rs0_table(i_loop,i_M_long_index)=Rs_cal(H_pf_total,H_pf_total,p);
        %% Calcu. rate with Imperfect channel with quan AoA but perfect grids
        [H_fb_quan] = For_Parfor_Calcu_ProposedScheme_with_perfect_AoD_Rate(L1,quan_Hss_rum_trans,A_total,B,N,M_long_max,K,spt_col_with_Mlong,UM_with_M_long)
        [~,D,THETA_quan]=ACE_RIS_2bit(N_iter,H_fb_quan,H_d);
        H_fb_quan_total= RIS_Change_Cascaded_to_Total_channel(H_fb_quan,H_d,THETA_quan);
        Hss_pf_quan_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_quan);
        %Cal Rate
        Rs1_table(i_loop,i_M_long_index)=Rs_cal(Hss_pf_quan_total,H_fb_quan_total,p);
        
        %% Imperfect channel with quan AoA and on-grids
        Value_M_long_table_of_i_M_long_index=M_long_table(i_M_long_index);        
   [H_fb_quan_on_gird] = For_Parfor_Calcu_ProposedScheme_with_Gt_resolution_AoD_Rate(L1,quan_on_gird_Hss_rum_trans,A_total,B,N,Value_M_long_table_of_i_M_long_index,K,spt_col_with_Mvar,UM_with_M_var);     
        [~,D,THETA_quan_on_gird]=ACE_RIS_2bit(N_iter,H_fb_quan_on_gird,H_d);
        H_fb_quan_on_gird_total= RIS_Change_Cascaded_to_Total_channel(H_fb_quan_on_gird,H_d,THETA_quan_on_gird);
        Hss_pf_quan_on_gird_total = RIS_Change_Cascaded_to_Total_channel(Hss,H_d,THETA_quan_on_gird);
        %Cal Rate
        Rs2_table(i_loop,i_M_long_index)=Rs_cal(Hss_pf_quan_on_gird_total,H_fb_quan_on_gird_total,p);
        
%         %
        
        
        fprintf('For FIG5:i_loop=%d of %d,i_M_long_index=%d\n',i_loop,loop_num,i_M_long_index);
        
    end
end

filename   =   strcat('Rs0_table',   '.mat');
save(['.\Data\FIG5data\',filename],    'Rs0_table')%,'-v7.3');

filename   =   strcat('Rs1_table',   '.mat');
save(['.\Data\FIG5data\',filename],    'Rs1_table')%,'-v7.3');

filename   =   strcat('Rs2_table',   '.mat');
save(['.\Data\FIG5data\',filename],    'Rs2_table')%,'-v7.3');
mRs0=mean(Rs0_table,1);
mRs1=mean(Rs1_table,1);
mRs2=mean(Rs2_table,1);
figure;
plot(M_long_index_table(1:length(M_long_index_table)-1),mRs0(1:length(M_long_index_table)-1),'k--','LineWidth',1.5);
hold on;
hold on;
plot(M_long_index_table(1:length(M_long_index_table)-1),mRs1(1:length(M_long_index_table)-1),'rd-','LineWidth',1.5);
hold on;
plot(M_long_index_table(1:length(M_long_index_table)-1),mRs2(1:length(M_long_index_table)-1),'bd-','LineWidth',1.5);
hold on;

set(gca,'XTick',M_long_index_table)
set(gca,'xticklabel',{'16','32','64','128','256','512','1024'});
legend('Perfect CSIT','Proposed scheme with perfect AoDs ','Proposed scheme with imperfect AoDs');
% axis([3 15 0 4]);
grid on;
xlabel('AoD resolution G_t at the BS ');
ylabel('Achievable sum-rate [bps/Hz]');





