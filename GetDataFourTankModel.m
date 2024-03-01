function [sys, constr, HankelMatrices] = GetDataFourTankModel(T, N, Tini, wvar, Ts)
    % GetDataFourTankModel returns all the required system information to
    % construct an MPC/DPC controller for it. Useful to swap systems around
    % in simulations without having to change much. 
    %
    % This is a linearized Four Tank model, 4th order discrete time TITO system following
    %   Berberich et al, 2021, "Data-driven model predictive control:
    %   closed-loop guarantees and experimental results"
    %
    % Note that this system is closed-loop unstable without any stability
    % criterion for prediction horizons less than N=70.
    %
    % arguments:
    %   T -- length of Hankel matrices
    %   N -- Prediction horizon
    %   Tini -- past window length
    %   wvar -- desired noise variance
    %   Ts -- Sampling period
    %
    % returns:
    %   sys -- ss object of system
    %   constr -- struct with the system constraints
    %   HankelMatrices -- struct with the Hankel matrices, also contains
    %                     original data

    A = [0.921 0 0.041 0;...
          0 0.918 0 0.033;...
          0 0 0.924 0;...
          0 0 0 0.937];
    B = [0.017 0.001; 0.001 0.023; 0 0.061; 0.072 0];

    C = [1 0 0 0; 0 1 0 0];
    D = [];

    sys = ss(A, B, C, D, Ts);
    n = size(A,1);
    nu = size(B,2);
    ny = size(C,1);

    umax = [2;2];
    umin = [-2;-2];
    ymax = [2;2];
    ymin = [-2;-2]; 

    constr = struct('umax', umax, 'umin', umin, 'ymax', ymax, 'ymin', ymin);

    %% generate data
    U = 1*idinput([T+N+Tini, nu], 'PRBS', [0, 1], [-1, 1])';
    X = zeros(n, size(U,2));
    Y = zeros(ny,size(U,2));

    %sensor noise
    W = wvar*randn(ny,size(U,2));

    for k=1:size(U,2)
        Y(:,k) = C*X(:,k) + W(:,k);
        if(k < size(U,2))
            X(:,k+1) = A*X(:,k) + B*U(:,k);
        end
    end
    disp(['SNR output noise: ', num2str(snr(Y, W))]);

    %% Define Controller
    Up = zeros(Tini*nu, T);
    Uf = zeros(N*nu, T);
    Yp = zeros(Tini*ny, T);
    Yf = zeros(N*ny, T);

    for i = 1:Tini
        Up((i-1)*nu+1:i*nu, :) = U(:, i  :i+T-1);
        Yp((i-1)*ny+1:i*ny, :) = Y(:, i+1:i+T);
    end

    for i = 1:N
        Uf((i-1)*nu+1:i*nu, :) = U(:, i+Tini  :i+Tini+T-1);
        Yf((i-1)*ny+1:i*ny, :) = Y(:, i+Tini+1:i+Tini+T);
    end

    HankelMatrices = struct('Up', Up, 'Uf', Uf, 'Yp', Yp, 'Yf', Yf, 'U', U, 'Y', Y);
end

