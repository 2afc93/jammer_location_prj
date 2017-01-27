% function [ x_state, P_cov_new, K_gain ] = EKFilter_new( x_uav_0, x_uav_k, alpha, x_jam_old, P_cov_old, F_KF, G_KF, Q_KF, R_KF )
function [ system_posteriori ] = EKFilter_new(system_priori, F_KF, G_KF, Q_KF, R_KF)
    % Extended Kalman Filter
    k = system_priori.step;
    x_jam_old = system_priori.x_jam(:,k-1);
    if k > 100
        R_KF = 0.07; 
        alpha = system_priori.alpha(k)/system_priori.alpha(k-5);
        x_uav_0 = system_priori.x_uav(:,k-5);
    else
        alpha = system_priori.alpha(k);
        x_uav_0 = system_priori.x_uav(:,1);
    end
    x_uav_k = system_priori.x_uav(:,k);
    P_cov_old = system_priori.P_cov_state(:,:,k-1);
    system_posteriori = system_priori;
    
    %% Prediction
   
    % Equation 1
    X_s = F_KF*x_jam_old;

    % Equation 2
    P_s = F_KF*P_cov_old*F_KF' + G_KF*Q_KF*G_KF'; % Error covariance propagation

    %% Correction

    % nonlinear measurement eq
    h = hk(x_uav_0, x_uav_k, X_s);

    % Jacobian of nonlinear measurement eq.
    H = JH(x_uav_0, x_uav_k, X_s);

    % Equation 3: Innovation
    v = (alpha - h);

    % % IMPORTANT!! - ANGLE MANIPULATION ==============
    % % make the absolute value of angles under 180deg
    % if abs(v(2))>pi
    %     v(2) = v(2) - 2*pi*sign(v(2));
    % end

    % Equation 4: Innovation matrix
    S = (H*P_s*H' + R_KF);

    % Equation 5: Kalman gain
    K_gain = P_s*H'/S;

    % Equation 6: State update
    x_state = X_s + K_gain*v;

    % Equation 7: Error covariance update
    P_cov_new = (eye(2) - K_gain*H)*P_s;
    
    system_posteriori.P_cov_state(:,:,k) = P_cov_new;
    system_posteriori.x_jam(:,k) = x_state;
    system_posteriori.K_gain(:,k) = K_gain;
    

    N_back = 450;
    alpha_filt = 0.7;
    if k < 1600
        system_posteriori.x_jam(:,k) = x_state;
    else
        system_posteriori.x_jam(:,k) = ( alpha_filt*x_state + (1-alpha_filt)*mean(system_posteriori.x_jam(:,(k-N_back):(k-1)),2) );
    end

end



% h(X): Nonlinear measurement eq
function h = hk(x_uav_0, x_uav_k, x_jam_k)
    h = norm([x_jam_k;0] - x_uav_0)^2/norm([x_jam_k;0] - x_uav_k)^2;
end

function h = h_alpha_ext( x_uav_0, x_uav_k, x_jam_k, Gt_ratio)
    h = Gt_ratio*norm([x_jam_k;0] - x_uav_0)^2/norm([x_jam_k;0] - x_uav_k)^2;
end

%% partial_h/partial_X: Jacobian of measurement eq
function H = JH(x_uav_0, x_uav_k, x_jam_k)
    H(1,1) = 2*(x_jam_k(1) - x_uav_0(1))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(1) - x_uav_k(1))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
    H(1,2) = 2*(x_jam_k(2) - x_uav_0(2))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(2) - x_uav_k(2))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
end

function H = JH_ext(x_uav_0, x_uav_k, x_jam_k, Gt_ratio)
    H(1,1) = Gt_ratio*(2*(x_jam_k(1) - x_uav_0(1))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(1) - x_uav_k(1))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4));
    H(1,2) = Gt_ratio*(2*(x_jam_k(2) - x_uav_0(2))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(2) - x_uav_k(2))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4));
    H(1,3) = norm([x_jam_k;0] - x_uav_0)^2/norm([x_jam_k;0] - x_uav_k)^2;
end

