function [ x_state, P_cov_new, K_gain ] = EKFilter( x_uav_0, x_uav_k, alpha, x_jam_old, P_cov_old, F_KF, G_KF, Q_KF, R_KF )

    % Extended Kalman Filter

    %% Prediction


    % Equation 1
    X_s = F_KF*x_jam_old;

    % Equation 2
    P_s = F_KF*P_cov_old*F_KF' + G_KF*Q_KF*G_KF'; % Error covariance extrapolation

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


end



% h(X): Nonlinear measurement eq
function h = hk(x_uav_0, x_uav_k, x_jam_k)
    h = norm([x_jam_k;0] - x_uav_0)^2/norm([x_jam_k;0] - x_uav_k)^2;
end


%% partial_h/partial_X: Jacobian of measurement eq
function H = JH(x_uav_0, x_uav_k, x_jam_k)
    H(1,1) = 2*(x_jam_k(1) - x_uav_0(1))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(1) - x_uav_k(1))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
    H(1,2) = 2*(x_jam_k(2) - x_uav_0(2))/(norm([x_jam_k;0] - x_uav_k)^2) - (2*(x_jam_k(2) - x_uav_k(2))*norm([x_jam_k;0] - x_uav_0)^2)/(norm([x_jam_k;0] - x_uav_k)^4);
end

