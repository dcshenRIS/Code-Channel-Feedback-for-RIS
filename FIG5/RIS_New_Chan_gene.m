function [Hss,quan_Hss_rum_trans,quan_on_gird_Hss_rum_trans,A_quan,spt_col_with_Mlong,spt_col_with_Mvar]=RIS_New_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_var,M_long)
%% descriptions of this function
% This function is to gene. the cascaded channel. Note that the difference
% between this fun. and RIS_New_offgird_Chan_gene.m is that, this fun. is
% used to gene. the channel with on-grid AoD and continuous AoAs. While the
% latter one is used to gene. the channel with continuous AoD and continuous AoAs.
% Specifically, to make the simulation more convenient, we gene. the perfect AoDs by setting AoD resolution as
% M_long_max.The perfect AoDs (the red curve) are selected randomly from
% the 2048 grids.  Then, limited by the quantization accuracy, these AoDs
% are quantized with a lower resolution from 16 to 1024.
% ---------------- input descriptions -------------------------------------
%See in RIS_New_offgird_Chan_gene.m
% ---------------- output descriptions ------------------------------------
%See in RIS_New_offgird_Chan_gene.m


P=L2;
d=0.5;

G=zeros(M,N);
H=zeros(N,N,K);
ALPHA1=zeros(L2,K);%channel path gain for RIS-UE channel 
ALPHA2=zeros(L1,1);%channel path gain for BS-RIS channel 
A1=zeros(N,K,L2);%Steering matrix for 'AoD' at the RIS
A2=zeros(N,L1);%Steering matrix for 'AoA' at the RIS
A1_quan=zeros(N,K,L2);%每个用户RIS离开角的导引矢量(基于量化角度生成)
A2_quan=zeros(N,L1);%RIS到达角的导引矢量(基于量化角度生成)

%% BS AoD格点化
Mindex=[-(M_long-1)/2:1:(M_long/2)]'*(2/M_long);
rand_index_Mindex=randperm(M_long);
rand_Mindex=Mindex(rand_index_Mindex);

spt_col_with_Mlong=M_long+1-rand_index_Mindex(1:L1);
psi=rand_Mindex(1:L1);

Mindex_m_var=[-(M_var-1)/2:1:(M_var/2)]'*(2/M_var);%
[trans_on_gird_psi_with_Mvar,on_gird_psi_index_with_Mvar] = On_gird2Off_gird(psi,Mindex_m_var);
psi_on_gird_with_Mvar=trans_on_gird_psi_with_Mvar';
spt_col_with_Mvar=M_var+1-on_gird_psi_index_with_Mvar(1:L1);


%% RIS到达角非格点化
phi21=random('unif',-1,1,1,L1);
phi21_quan=UQ(phi21,B0,-1,1);


phi22=random('unif',-1,1,1,L1);
phi22_quan=UQ(phi22,B0,-1,1);

%% generate G channel M*N
alpha2 = zeros(L1,1);
alpha2(1:L1) = rand(L1,1)+1j*rand(L1,1);
alpha2(1:L1) = alpha2./sqrt(2);
while (find(abs(alpha2)<0.01))
    alpha2(1:L1) = rand(L1,1)+1j*rand(L1,1);
    alpha2(1:L1) = alpha2./sqrt(2);
end
ALPHA2(:,1)=alpha2;

for l = 1:L1
    a21 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi21(l));
    a22 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi22(l));
    a2=kron(a21,a22);%RIS到达
    A2(:,l)=a2;   %N*L1
       
    a21_quan = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi21_quan(l));
    a22_quan = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi22_quan(l));
    a2_quan=kron(a21_quan,a22_quan);%RIS到达
    A2_quan(:,l)=a2_quan;    
   
    b  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*psi(l));
    B_steering(:,l)=b; %M*L1 
    b_on_gird_with_M_var  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*psi_on_gird_with_Mvar(l));
    B_steering_on_gird_with_M_var(:,l)=b_on_gird_with_M_var; %M*L1 
    
    G = G + alpha2(l)*b*a2';% M*N
end
G=sqrt(M*N/L1)*G;

%% Channel between RIS & UE
%H_data=zeros(N,L1,K);
for k=1:K
    hr=zeros(N,1);
    alpha1 = zeros(L2,1);%L2 path:UE to RIS
    alpha1(1:L2) = rand(L2,1)+1j*rand(L2,1);
    alpha1(1:L2) = alpha1./sqrt(2);
    while (find(abs(alpha1)<0.01))
        alpha1(1:L2) = rand(L2,1)+1j*rand(L2,1);
        alpha1(1:L2) = alpha1./sqrt(2);
    end
    
    ALPHA1(:,k)=alpha1;
    
    %% 非格点化角度
    phi11=random('unif',-1,1,1,L2);
    phi11_quan=UQ(phi11,B0,-1,1);
    %[on_gird_phi11,on_gird_phi11_index] = On_gird2Off_gird(phi11,N1index);
    phi12=random('unif',-1,1,1,L2);
    phi12_quan=UQ(phi12,B0,-1,1);
    %[on_gird_phi12,on_gird_phi12_index] = On_gird2Off_gird(phi12,N2index);
    for l = 1:L2
        a11 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi11(l));
        a12 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi12(l));
        a1=kron(a11,a12);
        A1(:,k,l)=a1;
        a11_quan = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi11_quan(l));
        a12_quan = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi12_quan(l));
        a1_quan=kron(a11_quan,a12_quan);
        A1_quan(:,k,l)=a1_quan;        
    %     a1 = 1/sqrt(N)*exp(-1i*2*pi*[0:N-1]'*d*phi11(l));
        hr = hr + alpha1(l)*a1;
        
    end
    hr=sqrt(N/L2)*hr;
%     hr=(normrnd(0, 1, N, 1) + 1i*normrnd(0, 1, N, 1)) / sqrt(2);
    H(:,:,k)=diag(hr');
end

Hss=zeros(N,M,K);
Hss_quan=zeros(N,M,K);
Hss_quan_on_grid=zeros(N,M,K);
on_gird_Hss=zeros(N,M,K);

quan_on_gird_Hss_rum=zeros(N,M_var,K);
UM_with_M_var=exp(1i*2*pi*[0:M-1]'*d*[-(M_var-1)/2:1:(M_var/2)]*(2/M_var));       %DFT矩阵 M*M
UM_with_M_var=(1/sqrt(M))*UM_with_M_var;

quan_Hss_rum=zeros(N,M_long,K);
UM_with_M_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFT矩阵 M*M_long
UM_with_M_long=(1/sqrt(M))*UM_with_M_long;

Haa=zeros(N,M,K);
a=0;
for i_K=1:1:K
    for i_L1=1:1:L1
        for i_L2=1:1:L2
            a=a+1;
            %N*M*K
            Hss(:,:,i_K)=Hss(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1(:,i_K,i_L2))')*A2(:,i_L1)*(B_steering(:,i_L1))';
            Hss_quan(:,:,i_K)=Hss_quan(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1)*(B_steering(:,i_L1))';
            Hss_quan_on_grid(:,:,i_K)=Hss_quan_on_grid(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1)*(B_steering_on_gird_with_M_var(:,i_L1))';
            
            quan_Hss_rum(:,spt_col_with_Mlong(i_L1),i_K)=quan_Hss_rum(:,spt_col_with_Mlong(i_L1),i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1);
            quan_on_gird_Hss_rum(:,spt_col_with_Mvar(i_L1),i_K)=quan_on_gird_Hss_rum(:,spt_col_with_Mvar(i_L1),i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1);
            
            A_quan(:,i_L2,i_L1,i_K)=ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1);
           
        end%
    end
    NMSE(i_K)=10*log10(norm(Hss_quan(:,:,i_K)-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2);
end


quan_Hss_rum_trans=quan_Hss_rum(:,spt_col_with_Mlong,:);
quan_on_gird_Hss_rum_trans=quan_on_gird_Hss_rum(:,spt_col_with_Mvar,:);

end
