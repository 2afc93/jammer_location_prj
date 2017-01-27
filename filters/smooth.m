%% Problem 9.3.4.2.1 Fixed-Lag Smoothing (Example 9.19)
% Z. Tian
% 5/7/14

clear all;
clc;
close all;
% AR(2) model parameter:

p=2;
q=0;
N=1;% #of observations
L_vector=0:10;

% representtive values of (a1,a2):
A1=[-1.5 -1 -0.75 0.5];
A2=[0.6 0.7 0.5 0.8];
b0=1;

SNR_dB=10;
SNR=10^(SNR_dB/10);
sigma_w2=1;
sigma_s2=SNR*sigma_w2;

K=50;% #of Kalman iterations
P11=zeros(length(L_vector),length(A1));

for ind_a=1:length(A1)
    a1 = A1(ind_a);
    a2 = A2(ind_a);
    
    %check condition for being a realizable minimum phase system: eqn (9.45)
    if a2>-1 && a2<1 && a2-a1>-1 && a2+a1>-1
        disp('chosen (a_1,a_2) is fine');
    else
        disp('error: condition for being minimum phase system eqn is violated!');
        beep
        beep
        beep
        pause;
    end
    
    c1 = (1+a2)/(1+a2-a1^2-a2^2+a1^2*a2-a2^3);
    sigma_u2 = SNR*sigma_w2/c1;
    
    for ind_l=1:length(L_vector)
        L=L_vector(ind_l);
        
        % state space model parameters
        F=[-a1 -a2;1 0];
        F_eig_mag=abs(eig(F));
        if F_eig_mag>1
            disp('Warning: the selected F matrix has eigen value greater than 1');
        end
        
        G=[b0;0];
        C=[1 0];
        Q=sigma_u2;
        R=sigma_w2;
        
        %% augmented state model parameters:
        F_aug=[F zeros(p,(L)*p);eye((L)*p) zeros((L)*p,p)];
        G_aug=[G;zeros((L)*p,1)];
        C_aug=[C zeros(1,(L)*p)];
        Q_aug=sigma_u2;
        R_aug=sigma_w2;
        
        p1=(L+1)*p;
        
        % step1: initial conditions:
        P0=sigma_s2*eye(p1);
        x0=zeros(p1,1);
        
        % Kalman estimates:
        x_p=zeros(p1,1,K);
        P_p=zeros(p1,p1,K);
        r_p=zeros(N,1,K);
        r_tilde=zeros(N,1,K);
        P_tilde=zeros(N,N,K);
        K_k=zeros(p1,N,K);
        x=zeros(p1,1,K);
        P=zeros(p1,p1,K);

        x(:,:,1)=x0;
        P(:,:,1)=P0;

        % actual process:
        x_actual=zeros(p1,1,K);
        %observation:
        r=zeros(N,1,K);
        
        for k=2:K
            % process:
            u=sqrt(sigma_u2)*randn(1,1);
            x_actual(:,:,k)=F_aug*x_actual(:,:,k-1)+G_aug*u;
            
            % observation:
            w=sqrt(sigma_w2)*randn(1,1);
            r(:,:,k)=C_aug*x_actual(:,:,k)+w;
        end
        
        for k=2:K
            % step2: prediction
            x_p(:,:,k)=F_aug*x(:,:,k-1);
            P_p(:,:,k)=F_aug*P(:,:,k-1)*F_aug'+G_aug*Q_aug*G_aug';
            
            r_p(:,:,k)=C_aug*x_p(:,:,k);
            r_tilde(:,:,k)=r(:,:,k)-r_p(:,:,k);
            P_tilde(:,:,k)=C_aug*P_p(:,:,k)*C_aug'+R_aug;
            
            % step3: Kalman gain
            K_k(:,:,k)=P_p(:,:,k)*C_aug'*inv(P_tilde(:,:,k));
            
            %step4: state and error covariance estimates
            x(:,:,k)=x_p(:,:,k)+K_k(:,:,k)*r_tilde(:,:,k);
            P(:,:,k)=(eye(p1)-K_k(:,:,k)*C_aug)*P_p(:,:,k);
            
        end
        P11(ind_l,ind_a)=squeeze(P(end,end,end))/sigma_s2;
        
    end
    
end


% plot results:
%% a) plot P vs L:
figure(1);
set(0,'DefaultLineLineWidth', 2, 'DefaultAxesFontSize',14, 'DefaultTextFontSize',14);
plot(L_vector,P11(:,1),':rd');
hold on;
plot(L_vector,P11(:,2),'-.m^');
hold on;
plot(L_vector,P11(:,3),'-bs');
hold on; 
plot(L_vector,P11(:,4),'--go');

xlabel('L');
ylabel('P_{L\infty}/\sigma^2_s');
title(['SNR=',num2str(SNR_dB),' dB']);
grid on;
legend('a_1=-1.5, a_2=0.6','a_1=-1, a_2=0.7','a_1=-0.75, a_2=0.5','a_1=0.5, a_2=0.8','Location', 'SouthEast');