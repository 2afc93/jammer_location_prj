function [ system_post ] = UKEK_smooth(system_prio, F, G, Q, R, N_back)
    % [system_post] = UKSmoother(system_prio, F, G, Q, R, N_back) returns a
    % state x and covariance P estimations on time k based on estimations made 
    % by forward UKF up to k and those done by the backward UKF starting at k+N_back. 
    
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
    
       
end

function [y_t_prio, S_t_prio] = smooth(system, t, F, Q, R)
    if system.step == t
        y_t_prio = 0;
        S_t_prio = 0;
    else
        
        [y_tp1_prio, S_tp1_prio] = smooth( system, t+1, F, Q, R );
        
        % Measurement Update
        z_tp1 = system.alpha(t+1);
        H_tp1 = JH(system.x_uav(:,1), system.x_uav(:,t+1), system.x_jam(:,t+1));
        y_tp1_post = y_tp1_prio + H_tp1'/R*z_tp1;
        S_tp1_post = S_tp1_prio + H_tp1'/R*H_tp1;
        
        % Time Update
        n = size(F);
        Q_inv = eye(n)/Q;
        K = S_tp1_post/(S_tp1_post + Q_inv);
        S_t_prio = F'*(eye(n) - K)*S_tp1_post*F;
        y_t_prio = F'*(eye(n) - K)*y_tp1_post;
    end
end


%% partial_h/partial_X: Jacobian of measurement eq
function H = JH(x_uav_0, x_uav_k, x_jam_k)
    H(1,1) = 2*(x_jam_k(1) - x_uav_0(1))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(1) - x_uav_k(1))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
    H(1,2) = 2*(x_jam_k(2) - x_uav_0(2))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(2) - x_uav_k(2))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
end

