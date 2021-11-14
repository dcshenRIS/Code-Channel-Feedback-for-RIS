function [S1,D,THETA]=ACE_RIS_2bit(N_iter,Joint_Refl_H,H_d)
%% descriptions of this function
% This function is for the CEO algorithm for jointly design the beamforming
% at the BS and the RIS. Note that the RIS is 2 bits phase shift.
% ---------------- input descriptions -------------------------------------
%   "N_iter" is the number of iter.
%   "Joint_Refl_H" is the cascaded channel 
%   "H_d" is the direct channel M*U
% 
% ---------------- output descriptions ------------------------------------
%   "THETA" is the phase shifting paremeters for RIS with 1*N.

%% Para.
N = size(Joint_Refl_H,1);%RIS单元数
M = size(Joint_Refl_H,2);%BS单元数
U = size(Joint_Refl_H,3);%用户数

rho = 0.2;
KK = 100;
% N = size(G,1);
% U = size(Hr,1);
p(1,:,:) = (1/4)*ones(1,N,3);
R = ceil((1-rho)*KK);
S = zeros(1,KK);
for iter = 1:N_iter
    for i = 1:N
        X(:,i) = randsrc(KK,1,[1,1j,-1,-1j; p(iter,i,1),p(iter,i,2),p(iter,i,3),1-sum(p(iter,i,:))]);
    end
    temp = 0;
    for k = 1:KK
        %A = diag(X(k,:));
        A = X(k,:);%1*N
        for i_u=1:1:U
            H_eff(i_u,:)=A*Joint_Refl_H(:,:,i_u);%U*M
        end
        H_eff=H_eff+H_d';

        %% Digital beamforming 
        D = H_eff'*inv(H_eff*H_eff');
        D = D./repmat(sqrt(sum(abs(D).^2,1)),M,1); %M*U
        
        %temp = norm(D,2);
        
        
        beta = sqrt(U/trace(D*D'));
        H_eq = H_eff*D;
        temp=0;
        SNR=1;
        for u=1:U
            sum_inf = sum(abs(H_eq(:,u)).^2)-abs(H_eq(u,u))^2;
            temp = temp + log2(1+abs(H_eq(u,u))^2/(sum_inf+U/(SNR*beta^2)));
        end
        %
        %         F = A*D;
        %         H_eq = H*F;
        %         beta = sqrt(KK/trace(F*F'));
        %         temp = 0;
        %         for u=1:U
        %             sum_inf = sum(abs(H_eq(:,u)).^2)-abs(H_eq(u,u))^2;
        %             temp = temp + log2(1+abs(H_eq(u,u))^2/(sum_inf+KK/(SNR*beta^2)));
        %         end
        S(k) = temp;
    end
    [S1,order] = sort(S,'ascend');
    r(iter) = S1(R);
    candidate = find(S>r(iter));
    den = size(candidate,2);
    if den == 0
        break;
    end
    
    for j = 1:N
        num1 = 0; num2=0; num3=0;%num4=0;
        for m = 1:den
            if X(candidate(m),j) == 1
                num1 = num1 + 1;
            elseif X(candidate(m),j) == 1j
                num2 = num2 + 1;
            elseif X(candidate(m),j) == -1
                num3 = num3 + 1;
            end
            %            num4 = den - num1 - num2 - num3;
            %            p_4 = 1 - p(iter,j,1) - p(iter,j,2) - p(iter,j,3);
        end
        p(iter+1,j,1) = num1/den;
        p(iter+1,j,2) = num2/den;
        p(iter+1,j,3) = num3/den;
        W = sum(p(iter+1,j,:));
        %         p(iter+1,j,1) = (num1*(p_4+p(iter,j,1)))/(num1+num4);
        %         p(iter+1,j,2) = (num2*(p_4+p(iter,j,2)))/(num2+num4);
        %         p(iter+1,j,3) = (num3*(p_4+p(iter,j,3)))/(num3+num4);
    end
    %     T = sum(S(candidate))/length(candidate);
    %     den = length(candidate);
    %     for j = 1:N
    %         num = 0;
    %         for m = 1:den
    %             if X(candidate(m),j) == 1
    %                 num = num + S(candidate(m))/T;
    %             end
    %         end
    %         p(iter+1,j) = num/den;
    %         if p(iter+1,j) > 1
    %             p(iter+1,j) = 1;
    %         end
    %     end
%fprintf('S=%f,i_iter=%d\n',S1(order(end)),iter); 
if iter==N_iter
    fprintf('boom!');
end
end
%fprintf('Finish=========================');


%D;
THETA = X(order(end),:);%1*N
% F = diag(X(order(end),:));
% H_ef=Hr*F*G+Hd;
% D=H_ef'*inv(H_ef*H_ef');
