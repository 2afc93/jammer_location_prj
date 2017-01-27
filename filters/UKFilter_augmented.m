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
    % Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
    % Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
    % pp.401-422, 2004. 
    %
    % By Yi Cao at Cranfield University, 04/01/2008
    %

function [ system_posteriori ] = UKFilter_augmented(system_priori, F_KF, G_KF, Q_KF, R_KF)
    %% Data Acquisition
    k                   = system_priori.step;                     
    x_state_old         = [system_priori.x_jam(:,k-1); system_priori.Gt_ratio(k-1)];
    alpha               = system_priori.alpha(k);
    x_uav_0             = system_priori.x_uav(:,1);
    x_uav_k             = system_priori.x_uav(:,k);
    P_cov_old           = system_priori.P_cov_state(:,:,k-1);
    system_posteriori   = system_priori;
       
    hmeas               = @(x_state_k)h_alpha_ext(x_state_k, x_uav_0, x_uav_k);     % Measurement function
    
    
    %% Prediction (Linear Process)
    % State Propagation
    x_pred = F_KF*x_state_old;

    % Equation 2
    P_pred = F_KF*P_cov_old*F_KF' + G_KF*Q_KF*G_KF';           % Error covariance extrapolation
    
    %% Measurement (Non Linear Measurement)
    
    % Sigma Points Selection
    L                   = numel(x_pred);                    % State space size
    m                   = numel(alpha);                     % Measurement Vector Size
    alpha_param         = 0.8;                              % default, tunable
    ki                  = 1;                                % default, tunable
    beta                = 2;                                % default, tunable
    lambda              = alpha_param^2*(L+ki)-L;           % scaling factor
    c                   = L+lambda;                         % scaling factor
    Wm                  = [lambda/c 0.5/c+zeros(1,2*L)];    % Mean Weights
    Wc                  = Wm;                               % Covariance Weights
    Wc(1)               = Wc(1)+(1-alpha_param^2+beta);     % weights for covariance
    c                   = sqrt(c);                          %
    
    X_sigmas = sigmas(x_pred, P_pred, c);
    [z_mean,Z_t,S,Pz] = ut(hmeas, X_sigmas, Wm, Wc, m, R_KF);   % Measurement Unscented Transformation
    Pxz = zeros(L,1);
    
    for jj = 1:2*length(x_pred)+1
        Pxz = Pxz + Wc(jj)*(X_sigmas(:, jj) - x_pred)*(Z_t(jj) - z_mean)';
    end
    
    % Innovation
    v = alpha - z_mean;                 
    K_gain = Pxz/S;                                         % Kalman Gain
    x_state = x_pred + K_gain*v;                            % State update
    P_cov_new = P_pred - K_gain*S*K_gain';                  % Covariance Update
   
    system_posteriori.P_cov_state(:,:,k) = P_cov_new;
    system_posteriori.x_jam(:,k) = x_state(1:2);
    system_posteriori.Gt_ratio(k) = x_state(3);
    system_posteriori.K_gain(:,k) = K_gain;
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

function [y,Y,P,Y1] = ut(h, X, Wm, Wc, n, R)
    %Unscented Transformation
    %Input:
    %        h: nonlinear transformation
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
    
    Y = h(X);
    y = sum(Y.*(ones(n,1)*Wm));
    
    Y1  = Y - y(:,ones(1,L));
    P   = Y1*diag(Wc)*Y1' + R;          
end

% h(X): Nonlinear measurement eq
function h = h_alpha(x_jam_k, x_uav_0, x_uav_k)
    h = ( (x_jam_k(1,:)-x_uav_0(1)).^2 + (x_jam_k(2,:)-x_uav_0(2)).^2 + x_uav_0(3).^2 ) ./( (x_jam_k(1,:)-x_uav_k(1)).^2 + (x_jam_k(2,:)-x_uav_k(2)).^2 + x_uav_k(3).^2 );
end

function h = h_alpha_ext(x_state_k, x_uav_0, x_uav_k)
    h = x_state_k(3,:).*(( (x_state_k(1,:)-x_uav_0(1)).^2 + (x_state_k(2,:)-x_uav_0(2)).^2 + x_uav_0(3).^2 ) ./( (x_state_k(1,:)-x_uav_k(1)).^2 + (x_state_k(2,:)-x_uav_k(2)).^2 + x_uav_k(3).^2 ));
end