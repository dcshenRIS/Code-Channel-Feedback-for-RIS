function [Hss,on_gird_Hss_rum_trans,spt_col,A,A_quan,g]=RIS_New_offgird_Chan_gene(M,N,N1,N2,L1,L2,K,B0,M_long)
%% descriptions of this function
% This function is to gene. cascaded channel.
% ---------------- input descriptions -------------------------------------
%   "M" is the BS antennas.
%   "N" is the UPA RIS with N1*N2.
%   "L1" is the number of BS-RIS channel paths
%   "L2" is the number of UE-RIS channel paths
%   "K" is the number of users
%   "B0" is the quan. bits for cascaded AoAs
%   "M_long" is the number of AoD resolutions at the BS
% ---------------- output descriptions ------------------------------------
%   "Hss" is the real channel with N*M*K.
%   "on_gird_Hss_rum_trans" is the hybrid channel (the mini version
%   tranmission in different matlab fun.) with N*L1*K. Note that
%   actually, hybrid channel is N*M_long*K, in which many columns are
%   zero-vectors. Hence, to speed the fun., we only transmit the L1
%   non-zero columns to the main fun.

%   "spt_col" is the index for non-zero column.
% Important:直接从导引矢量B_steering处记录下其对应的M中的索引，但由于B_steering在后面做了共轭转置，所以计算索引时从后往前算。
%   "A" is the steering matrix for cascaded AoA.
%   "A_quan" is the steering matrix for quan. cascaded AoA.
%   "g" is the channel path gain.





P=L2;
d=0.5;

%   NMSE=ones(1,K);
%  while(max(NMSE)>=1e-3)
%% Init
G=zeros(M,N);
H=zeros(N,N,K);

Ha=zeros(N,M,K);
Hs=zeros(N,M,K);
support_Col=zeros(L1,K);%每个用户非零元所在列的索引
support_Row=zeros(L1,L2,K);%每个用户非零元所在行的索引

A1=zeros(N,K,L2);%每个用户RIS离开角的导引矢量
A2=zeros(N,L1);%RIS到达角的导引矢量

A1_quan=zeros(N,K,L2);%每个用户RIS离开角的导引矢量(基于量化角度生成)
A2_quan=zeros(N,L1);%RIS到达角的导引矢量(基于量化角度生成)

%用来记录路径增益，以生成每个用户的路径增益向量（多用户构成矩阵）
ALPHA1=zeros(L2,K);%RIS-UE路径增益
ALPHA2=zeros(L1,1);%BS-RIS路径增益
GAIN=zeros(L1,L2,K);%从L1中给定一条BS-RIS路径，每个用户有L2条RIS-UE，因此有L2维路径增益。另外，G已经用了，故此处用GAIN。

%% 空间角度格点化
Mindex=[-(M-1)/2:1:(M/2)]'*(2/M);

Nindex=[-(N-1)/2:1:(N/2)]'*(2/N);
N1index=[-(N1-1)/2:1:(N1/2)]'*(2/N1);
N2index=[-(N2-1)/2:1:(N2/2)]'*(2/N2);

 %% DFT矩阵
% UM=exp(1i*2*pi*[0:M-1]'*d*[-(M-1)/2:1:(M/2)]*(2/M));       %DFT矩阵 M*M
% UM=(1/sqrt(M))*UM;

UN1=exp(1i*2*pi*[0:N1-1]'*d*[-(N1-1)/2:1:(N1/2)]*(2/N1));
UN2=exp(1i*2*pi*[0:N2-1]'*d*[-(N2-1)/2:1:(N2/2)]*(2/N2));
UN=kron(UN1,UN2);                        %DFT矩阵N*N
UN=(1/sqrt(N))*UN; 

UM_long=exp(1i*2*pi*[0:M-1]'*d*[-(M_long-1)/2:1:(M_long/2)]*(2/M_long));       %DFT矩阵 M*M_long
UM_long=(1/sqrt(M))*UM_long;
UM=UM_long;

%% generate G channel M*N
alpha2 = zeros(L1,1);
alpha2(1:L1) = rand(L1,1)+1j*rand(L1,1);
alpha2(1:L1) = alpha2./sqrt(2);
while (find(abs(alpha2)<0.01))
    alpha2(1:L1) = rand(L1,1)+1j*rand(L1,1);
    alpha2(1:L1) = alpha2./sqrt(2);
end

ALPHA2(:,1)=alpha2;%记录BS-RIS路径增益 L1*1

%% 原来生成格点化角度
rand_index_N1index=randperm(N1);
rand_N1index=N1index(rand_index_N1index);
phi21=rand_N1index(1:L1);
%phi21_quan=UQ(phi21,B0,-1,1);

rand_index_N2index=randperm(N2);
rand_N2index=N2index(rand_index_N2index);
phi22=rand_N2index(1:L1);
%phi22_quan=UQ(phi22,B0,-1,1);


rand_index_Mindex=randperm(M);
rand_Mindex=Mindex(rand_index_Mindex);
%=====================spt_col=这里可能没用====================================
spt_col=M+1-rand_index_Mindex(1:L1);
psi=rand_Mindex(1:L1);%取乱序角度的前L1个作为BS的L1条路径离开角2
psi_quan=UQ(psi,B0,-1,1);
%% 生成非格点话角度

N2index=[-(N2-1)/2:1:(N2/2)]'*(2/N2);
phi21=random('unif',-1,1,1,L1);
phi21_quan=UQ(phi21,B0,-1,1);
%[on_gird_phi21,on_gird_phi21_index] = On_gird2Off_gird(phi21,N1index);

phi22=random('unif',-1,1,1,L1);
phi22_quan=UQ(phi22,B0,-1,1);
%[on_gird_phi22,on_gird_phi22_index] = On_gird2Off_gird(phi22,N2index);

psi=random('unif',-1,1,1,L1);

Mindex_long=[-(M_long-1)/2:1:(M_long/2)]'*(2/M_long);
[on_gird_psi,on_gird_psi_index] = On_gird2Off_gird(psi,Mindex_long);%这里是格点角度（DOA）
spt_col=M_long+1-on_gird_psi_index(1:L1);
on_gird_psi;
NMSEANGLE=norm(psi-on_gird_psi)^2/norm(psi)^2;


for l = 1:L1
    a21 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi21(l));
    a22 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi22(l));
    a2=kron(a21,a22);%RIS到达
    A2(:,l)=a2;   %N*L1
    
%    on_gird_a21 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*on_gird_phi21(l));
%     on_gird_a22 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*on_gird_phi22(l));
%     on_gird_a2=kron(on_gird_a21,on_gird_a22);%RIS到达
%     on_gird_A2(:,l)=on_gird_a2;   %N*L1
    
    
    a21_quan = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*phi21_quan(l));
    a22_quan = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*phi22_quan(l));
    a2_quan=kron(a21_quan,a22_quan);%RIS到达
    A2_quan(:,l)=a2_quan;    
    
    on_gird_b  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*on_gird_psi(l));
    on_gird_B_steering(:,l)=on_gird_b; %M*L1
    
    b  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*psi(l));
    B_steering(:,l)=b; %M*L1    
    b_000  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*on_gird_psi(l));
    nmseb=norm(b_000-b)^2/norm(b)^2;
    
    b_quan  = 1/sqrt(M)*exp(-1i*2*pi*[0:M-1]'*d*psi_quan(l));
    B_quan(:,l)=b_quan;
%     b=kron(b1,b2);
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
    
    rand_index_N1index=randperm(N1);
    rand_N1index=N1index(rand_index_N1index);
    phi11=rand_N1index(1:L2);
    %phi11_quan=UQ(phi11,B0,-1,1);
    
    rand_index_N2index=randperm(N2);
    rand_N2index=N2index(rand_index_N2index);
    phi12=rand_N2index(1:L2);
    %phi12_quan=UQ(phi12,B0,-1,1);
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
        
%         on_gird_a11 = 1/sqrt(N1)*exp(-1i*2*pi*[0:N1-1]'*d*on_gird_phi11(l));
%         on_gird_a12 = 1/sqrt(N2)*exp(-1i*2*pi*[0:N2-1]'*d*on_gird_phi12(l));
%         on_gird_a1=kron(on_gird_a11,on_gird_a12);
%         on_gird_A1(:,k,l)=on_gird_a1;
        
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
on_gird_Hss=zeros(N,M,K);
on_gird_Hss_rum=zeros(N,M_long,K);
temp0904=zeros(N,M_long,K);

Haa=zeros(N,M,K);
a=0;
for i_K=1:1:K
    for i_L1=1:1:L1
        for i_L2=1:1:L2
            a=a+1;
            %N*M*K
            Hss(:,:,i_K)=Hss(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1(:,i_K,i_L2))')*A2(:,i_L1)*(B_steering(:,i_L1))';
%           错了on_gird_Hss(:,:,i_K)=on_gird_Hss(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((on_gird_A1(:,i_K,i_L2))')*on_gird_A2(:,i_L1)*(on_gird_B_steering(:,i_L1))';
            %on_gird_Hss(:,:,i_K)=on_gird_Hss(:,:,i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1(:,i_K,i_L2))')*A2(:,i_L1)*(on_gird_B_steering(:,i_L1))';
            on_gird_Hss_rum(:,spt_col(i_L1),i_K)=on_gird_Hss_rum(:,spt_col(i_L1),i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1(:,i_K,i_L2))')*A2(:,i_L1);
%A(:,i_L2,i_L2,i_K)=diag((A1(:,i_K,i_L2))')*A2(:,i_L1);
            A(:,i_L2,i_L1,i_K)=ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1(:,i_K,i_L2))')*A2(:,i_L1);
            A_quan(:,i_L2,i_L1,i_K)=ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1);
            g(i_L2,i_L1,i_K)=sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1);
            temp0904(:,spt_col(i_L1),i_K)=temp0904(:,spt_col(i_L1),i_K)+sqrt(N/L2)*sqrt(M*N/L1)*ALPHA1(i_L2,i_K)*ALPHA2(i_L1,1)*diag((A1_quan(:,i_K,i_L2))')*A2_quan(:,i_L1);
         
        end%
    end
    NMSE(i_K)=10*log10(norm(temp0904(:,:,i_K)*UM'-Hss(:,:,i_K),'fro')^2/norm(Hss(:,:,i_K),'fro')^2);
end


 %end
  NMSE;
%on_gird_Hss_rum_trans=zeros(N,L1,K); 
on_gird_Hss_rum_trans=on_gird_Hss_rum(:,spt_col,:);

end
