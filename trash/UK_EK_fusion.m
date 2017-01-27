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

function [ system_posteriori ] = UK_EK_fusion(system_priori, F_KF, G_KF, Q_KF, R_KF)
end

function [ x_new, P_new ] = EKF_range(x_old, P_old, alpha, t_new, F, Q, R, G)
    
    F = eye(2);
    Q = diag([4, 0.01]);
    R = 3^2;
    G = 1;
    
    %% Prediction
    % State Prediction
    x_s = F*x_old;

    % Covariance Propagation
    P_s = F*P_cov_old*F' + G*Q*G'; % Error covariance extrapolation
    
    % Measurement prediction
    H = [x_s(2) 0];
    y = H*x_s;
    
    %% Update
    % Measurement Innovation
    v = alpha(t) - y;
    
    % Innovation Matrix
    S = (H*P_s*H' + R);
    
    % Kalman Gain
    K = P_s*H'/S;
    
    % State Update
    x_new = x_s + K*v;
    
    P_new = P_s - K*H*P_s;
end

function [x_t, P_t] = rec_smoother(t, t_end, alpha, F, Q, R)
    
    if t > t_end
        x_t = 1;
        P_t = 0.01;
    else
        [x_tp1, P_tp1] = rec_smoother( t+1, t_end, alpha, F, Q, R );
        [x_t, P_t] = EKF(x_tp1, P_tp1, alpha, t, F, Q, R);
    end
end
    
function something()
%% Step definitions
    t_start = 26;               % This is the simulation step at which the filtering starts
    t_end   = system_prio.step;
    
    %% Forward Filtering: UKF estimation at simulation step  
    % The standard UKF filter is employed to get the estimation at the
    % simulation time
    
    
    
    system_post = UKFilter(system_prio, F, G, Q, R);
    t_smooth = t_end - N_back;
    if t_smooth > t_start
        [y_t_prio, S_t_prio] = smooth(system, t_smooth, F, Q, R);
        
        P_back_t = eye(size(S_t_prio))/S_t_prio;
        x_back_t = P_t_back*y_t_prio;
        
        x_for_t = system_post.x_jam(:,t_smooth);
        P_for_t = system_post.P_cov_state(:,:,t_smooth);
        
        K_t = P_for_t/(P_for_t + P_back_t);
        P_s_t = K_t*P_back_t;
        x_s_t = x_for_t + K_t*(x_back_t - x_for_t);
        
        system_post.x_jam(:,t_smooth) = x_s_t;
        system_post.P_cov_state(:,:,t_smooth) = P_s_t;
    end
    
    
    %% Prediction (Linear Process)
    % State Propagation
    x_pred = F_KF*x_state_old;

    % Equation 2
    P_pred = F_KF*P_cov_old*F_KF' + G_KF*Q_KF*G_KF';           % Error covariance extrapolation
    
    %% Measurement (Non Linear Measurement)
    
    % Sigma Points Selection
    L                   = numel(x_pred);                    % State space size
    m                   = numel(alpha);                     % Measurement Vector Size
    alpha_param         = 0.7;                              % default, tunable
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