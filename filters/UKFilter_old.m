% UKF   Unscented Kalman Filter for nonlinear dynamic systems
    % [x, P] = UKFilter(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
    % for nonlinear dynamic system (for simplicity, noises are assumed as additive):
    %           x_k+1 = f(x_k) + w_k
    %           z_k   = h(x_k) + v_k
    % where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
    %       v ~ N(0,R) meaning v is gaussian noise with covariance R
    % Inputs:   fstate: function handle for f(x)
    %           x: "a priori" state estimate
    %           P: "a priori" estimated state covariance
    %           h: fanction handle for h(x)
    %           z: current measurement
    %           Q: process noise covariance 
    %           R: measurement noise covariance

    % Output:   x: "a posteriori" state estimate
    %           P: "a posteriori" state covariance
    %           K: "a posteriori" Kalman Gain
    %
    %   [ x_state, P_cov_new, K_gain ] = UKFilter(x_uav_0, x_uav_k, alpha, x_jam_old, P_cov_old, F_KF, G_KF, Q_KF, R_KF) 
    %   returns the state estimate, x_state, state covariance, P_cov_new and Kalman Gain, K_gain
    
    % Inputs:   
    %           x_uav_0: first UAV position
    %           x_uav_k: last UAV position
    %           x_jam_old: "a priori" state estimate
    %           P_cov_old: "a priori" estimated state covariance
    %           h: fanction handle for h(x)
    %           alpha: current measurement
    %           Q_KF: process noise covariance 
    %           R_KF: measurement noise covariance
    
    % Example:
    %{
    n=3;      %number of state
    q=0.1;    %std of process 
    r=0.1;    %std of measurement
    Q=q^2*eye(n); % covariance of process
    R=r^2;        % covariance of measurement  
    f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
    h=@(x)x(1);                               % measurement equation
    s=[0;0;1];                                % initial state
    x=s+q*randn(3,1); %initial state          % initial state with noise
    P = eye(n);                               % initial state covraiance
    N=20;                                     % total dynamic steps
    xV = zeros(n,N);          %estmate        % allocate memory
    sV = zeros(n,N);          %actual
    zV = zeros(1,N);
    for k=1:N
      z = h(s) + r*randn;                     % measurments
      sV(:,k)= s;                             % save actual state
      zV(k)  = z;                             % save measurment
      [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
      xV(:,k) = x;                            % save estimate
      s = f(s) + q*randn(3,1);                % update process 
    end
    for k=1:3                                 % plot results
      subplot(3,1,k)
      plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
    end
    %}
    % Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
    % Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
    % pp.401-422, 2004. 
    %
    % By Yi Cao at Cranfield University, 04/01/2008
    %



% function [x, P, K] = UKFilter(fstate,x,P,hmeas,z,Q,R)
function [ x_state, P_cov_new, K_gain ] = UKFilter(x_uav_0, x_uav_k, alpha, x_jam_old, P_cov_old, F_KF, G_KF, Q_KF, R_KF)
    
    L=numel(x_jam_old);                                 %numer of states
    m=numel(alpha);                                 %numer of measurements
    alpha_param=1e-3;                                 %default, tunable
    ki=1;                                       %default, tunable
    beta=2;                                     %default, tunable
    lambda = alpha_param^2*(L+ki)-L;                    %scaling factor
    c = L+lambda;                                 %scaling factor
    Wm = [lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
    Wc = [lambda/c + (1 - alpha_param^2 + beta) 0.5/c+zeros(1,2*L)];
    c=sqrt(c);
    
    fstate  = @(x_jam)f_state(x_jam, F_KF, Q_KF, G_KF);
    hmeas   = @(x_jam)h_alpha(x_jam, x_uav_0, x_uav_k);
    
    X = sigmas(x_jam_old, P_cov_old, c);                                %sigma points around x
    [x1,X1,P1,X2] = ut(fstate, X, Wm, Wc, L, Q_KF);         %unscented transformation of process
    % X1=sigmas(x1,P1,c);                           %sigma points around x1
    % X2=X1-x1(:,ones(1,size(X1,2)));               %deviation of X1
    [z1,Z1,P2,Z2] = ut(hmeas,X1,Wm,Wc,m,R_KF);         %unscented transformation of measurments
    P12=X2*diag(Wc)*Z2';                            %transformed cross-covariance
%     K_gain = P12*inv(P2);                              % Kalman Gain
    K_gain = P12/P2;                              % Kalman Gain
    x_state = x1 + K_gain*(alpha - z1);                              %state update
    P_cov_new = P1 - K_gain*P12';                                %covariance update
end

function [y,Y,P,Y1] = ut(h, X, Wm, Wc, n, R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed sampling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);

for k=1:L                   
    Y(:,k) = h(X(:,k));       
    y = y + Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;          
end

function X = sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];
end

% h(X): Nonlinear measurement eq
function h = h_alpha(x_jam_k, x_uav_0, x_uav_k)
%     h = norm([x_jam_k;0] - x_uav_0)^2./norm([x_jam_k;0] - x_uav_k).^2;
    h = ( (x_jam_k(1)-x_uav_0(1)).^2 + (x_jam_k(2)-x_uav_0(2)).^2 + x_uav_0(3).^2 ) ./( (x_jam_k(1)-x_uav_k(1)).^2 + (x_jam_k(2)-x_uav_k(2)).^2 + x_uav_k(3).^2 );
    
end

function X_s = f_state(x_jam, F, Q, G)
    X_s = F*x_jam;% + G*Q*G';
end
