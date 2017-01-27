function [ system_posteriori ] = IEKFilter_augmented(system_priori, F_KF, G_KF, Q_KF, R_KF)
    tic
    % Extended Kalman Filter
    k = system_priori.step;
    x_state_old = [system_priori.x_jam(:,k-1); system_priori.Gt_ratio(k-1)];
    alpha = system_priori.alpha(k);
    x_uav_0 = system_priori.x_uav(:,1);
    x_uav_k = system_priori.x_uav(:,k);
    x_uav_km1 = system_priori.x_uav(:,k-1);
    P_cov_old = system_priori.P_cov_state(:,:,k-1);
    system_posteriori = system_priori;
    
    %% Prediction
    % Equation 1
    X_s = F_KF*x_state_old;

    % Equation 2
    P_s = F_KF*P_cov_old*F_KF' + G_KF*Q_KF*G_KF'; % Error covariance extrapolation
    
    h_fun = @(x_state)h_alpha_ext(x_uav_0, x_uav_k, x_state);
    H_fun = @(x_state)JH_ext(x_uav_0, x_uav_k, x_state);
    
    Xkn = X_s;
    diff = 1;
    tol = 0.01;
    iter = 0;
    while (diff > tol && iter < 5)
        iter = iter + 1;
        % nonlinear measurement eq
        h = h_fun(Xkn);
        % Jacobian of nonlinear measurement eq.
        H = H_fun(Xkn);
        
        % Equation 3: Innovation
        v = (alpha - h);
        % Equation 4: Innovation matrix
        S = (H*P_s*H' + R_KF);
        % Equation 5: Kalman gain
        K_gain = P_s*H'/S;
        
        Xkn_temp = X_s + K_gain*(v - H*(X_s - Xkn));
        
        diff = norm(abs(Xkn_temp - Xkn));
        
        Xkn = Xkn_temp;
    end
    diff
    % nonlinear measurement eq
    h = h_fun(Xkn);
    H = H_fun(Xkn);
    
    % Equation 3: Innovation
    v = (alpha - h);
    % Equation 4: Innovation matrix
    S = (H*P_s*H' + R_KF);
    % Equation 5: Kalman gain
    K_gain = P_s*H'/S;
    x_state = X_s + K_gain*(v - H*(X_s - Xkn));    
    % Error covariance update
    P_cov_new = (eye(3) - K_gain*H)*P_s;
    
    system_posteriori.P_cov_state(:,:,k) = P_cov_new;
    system_posteriori.x_jam(:,k) = x_state(1:2);
    system_posteriori.Gt_ratio(k) = x_state(3);
    system_posteriori.K_gain(:,k) = K_gain;
   
%     N_back = 450;
%     alpha_filt = 0.7;
%     if k < 1600
%         system_posteriori.x_jam(:,k) = x_state;
%     else
%         system_posteriori.x_jam(:,k) = ( alpha_filt*x_state + (1-alpha_filt)*mean(system_posteriori.x_jam(:,(k-N_back):(k-1)),2) );
%     end
%     system_posteriori.K_gain(:,k) = K_gain;
toc
end

function h = h_alpha_ext( x_uav_0, x_uav_k, x_state_k)
    h = x_state_k(3)*norm([x_state_k(1:2);0] - x_uav_0)^2/norm([x_state_k(1:2);0] - x_uav_k)^2;
end

function H = JH_ext(x_uav_0, x_uav_k, x_state_k)
    H(1,1) = x_state_k(3)*(2*(x_state_k(1) - x_uav_0(1))/(norm([x_state_k(1:2);0] - x_uav_k)^2) - (2*(x_state_k(1) - x_uav_k(1))*norm([x_state_k(1:2);0] - x_uav_0)^2)/(norm([x_state_k(1:2);0] - x_uav_k)^4));
    H(1,2) = x_state_k(3)*(2*(x_state_k(2) - x_uav_0(2))/(norm([x_state_k(1:2);0] - x_uav_k)^2) - (2*(x_state_k(2) - x_uav_k(2))*norm([x_state_k(1:2);0] - x_uav_0)^2)/(norm([x_state_k(1:2);0] - x_uav_k)^4));
    H(1,3) = norm([x_state_k(1:2);0] - x_uav_0)^2/norm([x_state_k(1:2);0] - x_uav_k)^2;
end