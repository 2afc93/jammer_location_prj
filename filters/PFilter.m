function [ system_posteriori ] = PFilter( system_priori, F_KF, G_KF, Q_KF, R_KF, Ns )
    
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
    global x_bnd y_bnd  % Global variables containing area dimensions
  
    k = system_priori.step;                 % Iteration number
    alpha = system_priori.alpha(k);         % Measurement at k-th iteration
    x_uav_0 = system_priori.x_uav(:,1);     % UAV initial position
    x_uav_k = system_priori.x_uav(:,k);     % UAV position at k-th iteration
    x_jam_old = system_priori.x_jam(:,k-1); % Jammer Previous Position
    particles = system_priori.particles;    % Particles population
    w         = system_priori.w;            % Particles Weights
    system_posteriori = system_priori;      % Output

    %% -----------------------------------------------------------------------------
    % FILTER INITIALISATION (k = 1)
    % ------------------------------------------------------------------------------
    % Particle sampling. Population defined as a random uniform distribution
    % within the search area.
    
%     if isempty(particles.x) || isempty(particles.y) || isempty(w)
    if any(isnan(particles.x(k-1,:))) || any(isnan(particles.y(k-1,:)))
        particles.x(k-1,:) = x_bnd/4 + x_bnd/2*rand(1,Ns);  % X coordinates
        particles.y(k-1,:) = y_bnd/4 + y_bnd/2*rand(1,Ns);  % Y coordinates
        w(:,k-1) = repmat(1/Ns, Ns, 1);                     % Equal Weights
    end
    
    % ------------------------------------------------------------------------------
    
    %% -----------------------------------------------------------------------------
    % PARTICLE FILTER (k > 1)
    % ------------------------------------------------------------------------------
    % Particle sampling at k-th iteration based on 1:k-1 previous
    % populations
    
    xkm1 = [particles.x(k-1,:); particles.y(k-1,:)]; % Previous Particle Population
    
    xk(1,:) = F_KF(1,:)*xkm1 + normrnd(0,sqrt(Q_KF(1,1)), 1, Ns);   % X sampling
    xk(2,:) = F_KF(2,:)*xkm1 + normrnd(0,sqrt(Q_KF(2,2)), 1, Ns);   % Y sampling
    
    % Weights computing
    wkm1 = w(:,k-1); % Previous Particles Weight
    v_hat = (alpha - h_alpha(xk, x_uav_0, x_uav_k))';   % Difference between computed and measured state
    
    if isnan(wkm1)
        wkm1 = repmat(1/Ns, Ns, 1);
        disp('Reassigning weights...')
    end
    % normpdf(X, MU, SIGMA) returns the gaussian distribution of mean MU
    % and standard deviation SIGMA evaluated at points in X.
    % wk = wkm1 .* (( 1 / sqrt(R_KF) / sqrt(2*pi) ) * exp( - v_hat.^2 / 2 / R_KF ));
    wk = wkm1 .* normpdf(v_hat, 0, sqrt(R_KF));     % Weights Vector
    wk = wk./sum(wk);                               % Weights Vector Normalisation

     % Resampling
    % Calculate effective sample size: eq 48, Ref 1
    Neff = 1/sum(wk.^2);                    % Efective Sample Size
    
    resample_percentage = 0.5;
    Nt = resample_percentage*Ns;
    
    if Neff < Nt
        [xk, wk, ~] = resample(xk, wk);
    end
    
    % Estimated state: the mean of the probability distribution
    system_posteriori.x_jam(:,k)        = xk*wk;

    % Save particle population and weights in k-th iteration
    system_posteriori.particles.x(k,:)  = xk(1,:); 
    system_posteriori.particles.y(k,:)  = xk(2,:);
    system_posteriori.w(:,k)            = wk;
    
    
    
end

% h(X): Nonlinear measurement eq
function h = h_alpha(x_jam_k, x_uav_0, x_uav_k)
    h = ( (x_jam_k(1,:)-x_uav_0(1)).^2 + (x_jam_k(2,:)-x_uav_0(2)).^2 + x_uav_0(3).^2 ) ./( (x_jam_k(1,:)-x_uav_k(1)).^2 + (x_jam_k(2,:)-x_uav_k(2)).^2 + x_uav_k(3).^2 );
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
