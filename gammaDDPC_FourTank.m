%% Intro
% Predictive control of a FourTank model using the
% gamma-- Data-driven Predictive Control (gDDPC) framework.
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

[Qs, Ls] = qr([Up;Yp;Uf;Yf]', 0); 
Ls = Ls'; Qs = Qs';
L11 = Ls(1:Tini*(ny+nu),                   1:Tini*(ny+nu));
L21 = Ls(Tini*(ny+nu)+1:Tini*(ny+nu)+N*nu, 1:Tini*(ny+nu));
L22 = Ls(Tini*(ny+nu)+1:Tini*(ny+nu)+N*nu, Tini*(ny+nu)+1:Tini*(ny+nu)+N*nu);
L31 = Ls(Tini*(ny+nu)+N*nu+1:end,          1:Tini*(ny+nu));
L32 = Ls(Tini*(ny+nu)+N*nu+1:end,          Tini*(ny+nu)+1:Tini*(ny+nu)+N*nu);
L33 = Ls(Tini*(ny+nu)+N*nu+1:end,          Tini*(ny+nu)+N*nu+1:Tini*(ny+nu)+N*nu+N*ny);

Q1 = Qs(1:Tini*(ny+nu),:);
Q2 = Qs(Tini*(ny+nu)+1:Tini*(ny+nu)+N*nu,:);
Q3 = Qs(Tini*(ny+nu)+N*nu+1:end,:);
L11i = pinv(L11);


%% Build controller
Q = 15*eye(ny);
R = 1e-1*eye(nu);
beta2 = 0;
beta3 = 30;

Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

% to plot regularization score, break below and inspect. Then re-run
%on the peak gave the best response!
% Ry = kron(ones(N,1), [0.65;0.77]);
% Ru = kron(ones(N,1), [0;0]);
% HankeRausRegularization([chol(Omega)*L33; chol(Psi)*zeros(N*nu,N*ny)], [chol(Omega)*Ry; chol(Psi)*Ru], [1e-4, 1e8, 1000]);

cpu=[];
u = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);
gamma2 = sdpvar(nu*N,1);
gamma3 = sdpvar(ny*N,1);

z_ini = sdpvar(Tini*(nu+ny), 1);

% Define objective function
objective = (y-ref)'*Omega*(y-ref)+(u)'*Psi*(u); % base MPC cost
objective = objective + beta2*(gamma2'*gamma2) + beta3*(gamma3'*gamma3);

% Build constraints
constraints = [u == L21*L11i*z_ini+L22*gamma2, y == L31*L11i*z_ini+L32*gamma2+L33*gamma3];

for k = 1:N   % and the bounds on the inputs and outputs
    constraints = [constraints, constr.umin<=u(nu*(k-1)+1:nu*k)<=constr.umax, constr.ymin<=y(ny*(k-1)+1:ny*k)<=constr.ymax];
end

Parameters = {z_ini, ref};
Outputs = {u, y};
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


rng(1904);
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
        Z_ini = [U_ini; Y_ini];

        [Sol, err] = controller({Z_ini, Rk});
        Uk = Sol{1}; 
        Yk = Sol{2}; %Yk can be used for evaluation, is not stored
        cpu = [cpu toc];
       
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
stairs(t, y(1,:), 'r', 'DisplayName', 'gDDPC y1');
stairs(t, y(2,:), 'b', 'DisplayName', 'gDDPC y2');
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
disp(['gDDPC ISE =', num2str(ISE), ' IAE = ', num2str(IAE), ' INEN = ', num2str(INEN), ' cpu = ', num2str(mean(cpu))]);

