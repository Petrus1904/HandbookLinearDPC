function [sys, constr, HankelMatrices] = GetDataPlaneModel(T, N, Tini, wvar, Ts)
    % GetDataPlaneModel returns all the required system information to
    % construct an MPC/DPC controller for it. Useful to swap systems around
    % in simulations without having to change much. 
    %
    % This system is a linearized plane model, 4th order TITO system in discrete time.
    % Model is presented in:
    % Camacho & Bordons, "Model Predictive Control", 2007
    %
    % The controller can control the plane's longitudinal velocity and climb
    % rate
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
    %   HankelMatrices -- struct with the Hankel matrices and the original
    %                     data sequences for sysID

    u0 = 7.74;
    Ac = [-0.003 0.039 0 -0.322; -0.065 -0.319 7.74 0; 0.02 -0.101 -0.429 0; 0 0 1 0];
    Bc = [0.010 1; -0.18 -0.04; -1.16 0.598; 0 0];

    Cc = [1 0 0 0; 0 -1 0 u0];
    D = [];

    boing_747c = ss(Ac,Bc,Cc,D,'StateName',{'Longitudinal Velocity(u)' 'Lateral Velocity(w)' 'Angular Velocity(w)' 'Angle(theta)'},...
                        'InputName',{'Elevator(e)' 'Throttle(t)'},...
                        'OutputName',{'Velocity(v)' 'Climb rate(h)'});

    sys = c2d(boing_747c,Ts,'zoh');
    A = sys.A;
    B = sys.B;
    C = sys.C;

    n = size(A,1);
    nu = size(B,2);
    ny = size(C,1);

    umax = [20;20];
    umin = [-20;-20];
    ymax = [25;15];
    ymin = [-25;-15]; 

    constr = struct('umax', umax, 'umin', umin, 'ymax', ymax, 'ymin', ymin);

    %% generate data
    U = 3*idinput([T+N+Tini, nu], 'PRBS', [0, 1], [-1, 1])';
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

