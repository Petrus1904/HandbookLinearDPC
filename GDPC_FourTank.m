%% Intro
% Predictive control of a four tank model using the
% Generalized Data-driven Predictive Control (GDPC) framework.
% Data acquisition and simulation are subject to sensor noise.

clear all;

%Set this to any constant to fix the noise seed every run.
rng("shuffle");

%% Define controller parameters and obtain Hankel matrices
Ts = 1;
N = 70;
Tini = 30;
T = 2500;

wvar = 0.02; %noise variance

[sys, constr, Hankel] = GetDataFourTankModel(T, N, Tini, wvar, Ts);

A = sys.A;
B = sys.B;
C = sys.C;

n = size(A,1);
nu = size(B,2);
ny = size(C,1);

Up = Hankel.Up;
Yp = Hankel.Yp;
Uf = Hankel.Uf;
Yf = Hankel.Yf;

Theta = Yf*pinv([Up;Yp;Uf]);
P1 = Theta(:,1:Tini*nu);
P2 = Theta(:, Tini*nu+1:Tini*(nu+ny));
Gamma = Theta(:, end-N*nu+1:end);

%% Build controller
Q =  15*eye(ny);
R =  1e-1*eye(nu);
l_g = 10;
l_y = 1e6;
Tg = 250;

Up_g = Up(:, 1:Tg);
Yp_g = Yp(:, 1:Tg);
Uf_g = Uf(:, 1:Tg);
Yf_g = Yf(:, 1:Tg);

Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);
Gi = inv(Psi + Gamma'*Omega*Gamma);

% to plot regularization score, break below and inspect. Then re-run
% Ry = kron(ones(N,1), [0.65;0.77]);
% Ru = zeros(N*nu, 1);
% uu = Gi*Gamma'*Omega*Ry;
% yy = Gamma*uu;
% HankeRausRegularization([chol(Omega)*Yf_g; chol(Psi)*Uf_g], [chol(Omega)*(Ry-yy); chol(Psi)*(Ru-uu)], [1e-4, 1e3, 1000], l_g);

cpu=[];
u = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);

g = sdpvar(Tg, 1);
u_spc = sdpvar(N*nu, 1);
y_spc = sdpvar(N*ny, 1);
sigma_y = sdpvar(Tini*ny, 1);

% Define objective function
objective = (y-ref)'*Omega*(y-ref)+(u)'*Psi*(u); % base MPC cost
objective = objective +l_g*(g'*g)+l_y*(sigma_y'*sigma_y); % cost for DeePC regularization terms

% Build constraints
constraints = [Up_g*g == zeros(Tini*nu,1), Yp_g*g == sigma_y, y-y_spc==Yf_g*g, u-u_spc==Uf_g*g];

for k = 1:N   % and the bounds on the inputs and outputs
    constraints = [constraints, constr.umin<=u(nu*(k-1)+1:nu*k)<=constr.umax, constr.ymin<=y(ny*(k-1)+1:ny*k)<=constr.ymax];
end

Parameters = {u_spc, y_spc, ref};
Outputs = {u, y, g};
options = sdpsettings('solver', 'osqp', 'verbose', 0, 'debug', 0);
controller = optimizer(constraints, objective, options, Parameters, Outputs);

%% Initialize Simulation
Tmax = 150;
t = 0:Ts:Tmax;
simLen = size(t,2);

x0 = zeros(n,1);
w = wvar*randn(ny, simLen);
%reference and disturbance sequences. Note that reference must be N samples
%longer to accommodate for the last predictions
r = [0.65; 0.77]*ones(1, simLen+N);
d = zeros(ny, simLen);


rng(1904); %this is only to fix the seed for the initial input window. The output noise is define before
y = zeros(ny, simLen);
u = zeros(nu, simLen);
x = zeros(n, simLen+1);
x(:,1) = x0;

%% Simulation
nbytes = fprintf('time: 0 of %d', Tmax);
err = 0;

for k = 1:simLen
    
    %Print current time and status without flooding the command window
    if(rem(k, 5)==0)
        fprintf(repmat('\b', 1, nbytes));
        nbytes = fprintf('processing %.3f of %d, QP status: %s', Ts*k, Tmax, yalmiperror(err));
    end
    
    %Compute new (measured) output
    y(:,k) = C*x(:,k) + w(:,k);

    if(k >= Tini+1)
        %Update control law
        Rk = r(:, k+1:k+N);
        Rk = Rk(:);
        tic;
        U_ini = u(:, k-Tini:k-1);
        U_ini = U_ini(:); %flatten vector
        Y_ini = y(:,k-Tini+1:k);
        Y_ini = Y_ini(:);
        
        Uspc = -Gi*Gamma'*Omega*(P1*U_ini + P2*Y_ini - Rk);
        Yspc = P1*U_ini + P2*Y_ini + Gamma*Uspc;

        [Sol, err] = controller({Uspc, Yspc, Rk});
        Uk = Sol{1}; 
        Yk = Sol{2}; %Yk can be used for evaluation, is not stored
        cpu = [cpu toc];
        if(err ~= 0)
            stop = 1; %add breakpoint here to evaluate what went wrong
        end
       
        u(:,k) = Uk(1:nu);
              
    else
        % To prevent initial infeasibility, let the system run in open loop
        % until enough data is gathered to build U_ini, Y_ini
       u(:,k) = 0.5*randn(nu, 1);
    end

    % Update system
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

%% Display simulation results

figure();
ax1 = subplot(211);
stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference 1');
hold on;
stairs(t, r(2,1:simLen), 'k', 'DisplayName', 'Reference 2');
stairs(t, y(1,:), 'r', 'DisplayName', ['GDPC y1, Tg=', num2str(Tg)]);
stairs(t, y(2,:), 'b', 'DisplayName', 'GDPC y2');
ylabel('Output [-]');
legend;
grid on; grid minor;
% axis([0 t(end) -1 15]);


ax3 = subplot(212);
hold on;
stairs(t, u(1,:), 'r');
stairs(t, u(2,:), 'b');
%stairs(t, u_ss(2)*ones(1,length(t)), '--k');
% axis([0 t(end) -20 20])
ylabel('Input [-]');
xlabel('Time [s]');
grid on; grid minor;

%%

CtrlCost = 0;
ISE = 0;
IAE = 0;
INEN = 0;
for k = Tini+1:simLen
    CtrlCost = CtrlCost + (y(:,k)-r(:,k))'*Q*(y(:,k)-r(:,k))+u(:,k)'*R*u(:,k);
    ISE = ISE + (y(:,k)-r(:,k))'*(y(:,k)-r(:,k));
    IAE = IAE + sum(abs(y(:,k)-r(:,k)));
    INEN = INEN + u(:,k)'*u(:,k);
end
disp(['GDPC ISE =', num2str(ISE), ' IAE = ', num2str(IAE), ' INEN = ', num2str(INEN), ' cpu = ', num2str(mean(cpu))]);
