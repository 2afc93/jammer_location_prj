function [ system_posteriori ] = PFilter_augmented( system_priori, F_KF, G_KF, Q_KF, R_KF, Ns )
    tic
    % --------------------------------------------------------------------------------
    % PFilter implements the basic particle filter algorithm with systematic 
    % resampling for the GPS Jammer Localization problem
    % --------------------------------------------------------------------------------
    % Programmed by:
    % Alvaro Fernandez (a.fernadez-cobo@cranfield.ac.uk)
    % Cranfield University, December 20, 2016
    % 
    % Adapted from: 
    % Diego Andres Alvarez Marin
    % Universidad Nacional de Colombia, February 29, 2012
    %
    % --------------------------------------------------------------------------------
    % Input parameters
    % --------------------------------------------------------------------------------
    %
    % system_priori: structure with system data up to k-1 step
    %   .step [scalar]: iteration number
    %   .alpha [1xNk]: Power ratio measurement. Nk is the total number of
    %   iterations
    %   .x_uav [2 x Nk]: UAV x and y postition coordinates at iteration k
    %   .particles: substructure containing the coordinates of the filter particles
    %       .x [Nk x Ns]: x coordinates of the Ns particles at the Nk
    %       iterations.
    %       .y [Nk x Ns]: y coordinates of the Ns particles at the Nk
    %       iterations.
    %   .w [Ns x Nk]: Weights of the Ns particles at the Nk iterations.
    % F_KF [2 x 2]: system state matrix. x(k) = F_KF*x(k-1)
    % G_KF [eye(1)]: noise matrix. (useless)
    % Q_KF [2 x 2]: Process Noise matrix
    % R_KF [scalar]: Measurement Noise
    % Ns [scalar]: Size of particle population
    %
    % --------------------------------------------------------------------------------
    % Output parameters
    % --------------------------------------------------------------------------------
    %
    % system_posteriori: same structure as system_priori
    %
    % --------------------------------------------------------------------------------
    % Reference:
    % [1] Arulampalam et. al. (2002).  A tutorial on particle filters for 
    %     online nonlinear/non-gaussian bayesian tracking. IEEE Transactions on 
    %     Signal Processing. 50 (2). p 174--188
    
    %% Working Variables definition 
   
    xy_bnd = 12000;                         % Search Area XY-Dimension
    
    
    k                   = system_priori.step;           % Iteration number
    alpha               = system_priori.alpha(k);       % Measurement at k-th iteration
    x_uav_0             = system_priori.x_uav(:,1);     % UAV initial position
    x_uav_k             = system_priori.x_uav(:,k);     % UAV position at k-th iteration
    particles           = system_priori.particles;      % Particles population
    wxy                 = system_priori.w.xy;           % Position Particles Weights
    wG                  = system_priori.w.G;            % Transmiter Gain Particle
    x_jam_old           = system_priori.x_jam(1:2,k-1); % Previous Jammer Position
    Gratio_old          = system_priori.Gt_ratio(k-1);  % Previous Gain transmitter ratio
    system_posteriori   = system_priori;                % Output

    %% -----------------------------------------------------------------------------
    % FILTER INITIALISATION (k = 1)
    % ------------------------------------------------------------------------------
    % Particle sampling. Population defined as a random uniform distribution
    % within the search area.
    
    min_xy = xy_bnd/4; max_xy = 3/4*xy_bnd;                     % Minimum/Maximum XY transmitter coordinates
    min_G = 0.85; max_G = 1.15;                            % Minimum/Maximum transmitter gain ratio
    
    if any(isnan(particles.x(k-1,:))) || any(isnan(particles.y(k-1,:))) || any(isnan(particles.G(k-1,:)))
        
        particles.x(k-1,:)  = min_xy + (max_xy - min_xy)*rand(1,Ns); % X coordinates particles
        particles.y(k-1,:)  = min_xy + (max_xy - min_xy)*rand(1,Ns); % Y coordinates particles
        wxy(:,k-1)          = repmat(1/Ns, Ns, 1);
        
        particles.G(k-1,:) = min_G + (max_G - min_G)*rand(1,Ns);    % Transmitter Gain ratio particles
        wG(:,k-1) = repmat(1/Ns, Ns, 1);                             % Equal Weights
    end
    
    % ------------------------------------------------------------------------------
    
    %% -----------------------------------------------------------------------------
    % PARTICLE FILTER (k > 1)
    % ------------------------------------------------------------------------------
    % Particle sampling at k-th iteration based on 1:k-1 previous
    % populations
    
    Gkm1 = particles.G(k-1,:);                                      % Previous Gain Particle Population
    wGkm1   = wG(:,k-1);
    if any(isnan(wGkm1)) || Gratio_old > 2
        Gkm1 = min_G + (max_G - min_G)*rand(1,Ns);                  % Reset Transmitter Gain ratio particles
        wGkm1 = repmat(1/Ns, Ns, 1);
        disp('Reassigning gain weights...')
    end
    
    Gk      = (F_KF(3,3)*Gkm1) + normrnd(0,sqrt(Q_KF(3,3)), 1, Ns);    % Gain ratio particles propagation
    
    v_hat   = (alpha - h_alpha_ext(repmat(x_jam_old,1,Ns), x_uav_0, x_uav_k, Gk))';
    
    wGk = wGkm1 .* normpdf(v_hat, 0, sqrt(R_KF));
    wGk = wGk/sum(wGk);
    
    Neff = 1/sum(wGk.^2);
    
    resample_percentage = 0.6;
    Nt = resample_percentage*Ns;
    
    if Neff < Nt
        [Gk, wGk, ~] = resample(Gk, wGk);
    end
    
    Gratio_new = Gk*wGk;
    
    
    xkm1 = [particles.x(k-1,:); particles.y(k-1,:)];                % Previous Position Particle Population
    
    % Position particle propagation
    xk(1,:) = F_KF(1,1:2)*xkm1 + normrnd(0,sqrt(Q_KF(1,1)), 1, Ns);   % X particles propagation
    xk(2,:) = F_KF(2,1:2)*xkm1 + normrnd(0,sqrt(Q_KF(2,2)), 1, Ns);   % Y particles propagation
    
    
    % Weights computing
    wxykm1 = wxy(:,k-1); % Previous Particles Weight
    v_hat = (alpha - h_alpha_ext(xk, x_uav_0, x_uav_k, repmat(Gratio_new,1,Ns)))';           % Difference between computed and measured state
    
    if isnan(wxykm1)
        wxykm1 = repmat(1/Ns, Ns, 1);
        disp('Reassigning weights...')
    end
    
    
    % normpdf(X, MU, SIGMA) returns the gaussian distribution of mean MU
    % and standard deviation SIGMA evaluated at points in X.
    % wk = wkm1 .* (( 1 / sqrt(R_KF) / sqrt(2*pi) ) * exp( - v_hat.^2 / 2 / R_KF ));
    wxyk = wxykm1 .* normpdf(v_hat, 0, sqrt(R_KF));         % Weights Vector
    wxyk = wxyk./sum(wxyk);                                 % Weights Vector Normalisation

     % Resampling
    % Calculate effective sample size: eq 48, Ref 1
    Neff = 1/sum(wxyk.^2);                            % Efective Sample Size
    
    resample_percentage = 0.6;
    Nt = resample_percentage*Ns;
    
    if Neff < Nt
        [xk, wxyk, ~] = resample(xk, wxyk);
    end    
    
    % Estimated state: the mean of the probability distribution
    system_posteriori.x_jam(:,k)        = xk*wxyk;

    % Save particle population and weights in k-th iteration
    system_posteriori.particles.x(k,:)  = xk(1,:); 
    system_posteriori.particles.y(k,:)  = xk(2,:);
    system_posteriori.particles.G(k,:)  = Gk;
    system_posteriori.w.xy(:,k)         = wxyk;
    system_posteriori.w.G(:,k)          = wGk;
    system_posteriori.Gt_ratio(k)       = Gratio_new;
    
    toc
    
end

% h(X): Nonlinear measurement eq
function h = h_alpha_ext( x_state, x_uav_k, x_uav_0, Gratio)
    % x_state = [x_jam_x; x_jam_y; Gt_ratio]
    h = Gratio.*(( (x_state(1,:)-x_uav_0(1)).^2 + (x_state(2,:)-x_uav_0(2)).^2 + x_uav_0(3).^2 ) ./( (x_state(1,:)-x_uav_k(1)).^2 + (x_state(2,:)-x_uav_k(2)).^2 + x_uav_k(3).^2 ));
end

% Resampling function
function [xk, wk, idx] = resample(xk, wk)
    % Systematic Resampling: latin hypercube sampling on wk
    Ns = length(wk);                % Ns = number of particles
    
    edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
    edges(end) = 1;                 % get the upper edge exact
    u1 = rand/Ns;
    
    % this works like the inverse of the empirical distribution and returns
    % the interval where the sample is to be found
%     [~, idx] = histc(u1:1/Ns:1, edges);
    [~,~,idx] = histcounts(u1:1/Ns:1, edges);
    xk = xk(:,idx);                 % extract new particles
    wk = 1/Ns*ones(Ns, 1);          % now all particles have the same weight
end
