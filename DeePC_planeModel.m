%% Intro
% Data-EnablEd Predictive Control of the plane model. Considers l2
% regularization on g.
% Data acquisition and simulation are subject to sensor noise.

clear all;

%Set this to any constant to fix the noise seed every run.
rng('shuffle');

%% Define controller parameters and obtain Hankel matrices
Ts = 0.1;
N = 20;
Tini = 20;
T = 250;

wvar = 0.25; %noise variance
%Put this to 2500 to give uniform results with the other controllers (noise
%wise)
[sys, constr, Hankel] = GetDataPlaneModel(2500, N, Tini, wvar, Ts);

A = sys.A;
B = sys.B;
C = sys.C;

n = size(A,1);
nu = size(B,2);
ny = size(C,1);

Up = Hankel.Up(:,1:T); %(:, 1:nu*(Tini+N):end);
Yp = Hankel.Yp(:,1:T); %(:, 1:ny*(Tini+N):end);
Uf = Hankel.Uf(:,1:T); %(:, 1:nu*(Tini+N):end);
Yf = Hankel.Yf(:,1:T); %(:, 1:ny*(Tini+N):end);

%% Build controller
Q =  10*eye(ny);
R =  0.01*eye(nu);
l_y = 1e8;
l_g = 1e3;

Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

% to plot regularization score, break below and inspect. Then re-run
% Wp0 = eye(T)-pinv([Up;Yp;Uf])*[Up;Yp;Uf];
% Ry = kron(ones(N,1), [10;0]);
% Ru = 0*kron(ones(N,1), u_ss);
% HankeRausRegularization([chol(Omega)*Yf; chol(Psi)*Uf]*Wp0, [chol(Omega)*Ry; chol(Psi)*Ru], [1e-2, 1e8, 1000]);
% HankeRausRegularization([chol(Omega)*Yf]*Wp0, [chol(Omega)*Ry], [1e-2, 1e8, 1000]);

cpu=[];
u = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);

g = sdpvar(T, 1);
sigma_y = sdpvar(ny*Tini,1);
u_ini = sdpvar(Tini*nu, 1);
y_ini = sdpvar(Tini*ny, 1);

% Define objective function
objective = (y-ref)'*Omega*(y-ref)+(u)'*Psi*(u); % base MPC cost
objective = objective +l_g*(g'*g) + l_y*(sigma_y'*sigma_y); %cost for DeePC regularization terms

% Build constraints
constraints = [u_ini==Up*g, y_ini+sigma_y==Yp*g, y==Yf*g, u==Uf*g];

for k = 1:N   % and the bounds on the inputs and outputs
    constraints = [constraints, constr.umin<=u(nu*(k-1)+1:nu*k)<=constr.umax, constr.ymin<=y(ny*(k-1)+1:ny*k)<=constr.ymax];
end

Parameters = {u_ini, y_ini, ref};
Outputs = {u, y, g};
options = sdpsettings('solver', 'osqp', 'verbose', 0, 'debug', 0); %, 'osqp.eps_abs', 1e-8, 'osqp.eps_rel', 1e-8);
controller = optimizer(constraints, objective, options, Parameters, Outputs);

%% Initialize Simulation
Tmax = 10;
t = 0:Ts:Tmax;
simLen = size(t,2);

x0 = zeros(n,1);
w = wvar*randn(ny, simLen);
%reference and disturbance sequences. Note that reference must be N samples
%longer to accommodate for the last predictions
r = [10; 0]*ones(1, simLen+N);
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
        U_ini = u(:, k-Tini:k-1);
        U_ini = U_ini(:); %flatten vector
        Y_ini = y(:,k-Tini+1:k);
        Y_ini = Y_ini(:);
 
        tic;

        [Sol, err] = controller({U_ini, Y_ini, Rk});
        Uk = Sol{1}; 
        Yk = Sol{2}; %Yk can be used for evaluation, is not stored
        cpu = [cpu toc];
       
        u(:,k) = Uk(1:nu);
              
    else
        % To prevent initial infeasibility, let the system run in open loop
        % until enough data is gathered to build U_ini, Y_ini
       u(:,k) = 1*randn(nu, 1);
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
stairs(t, y(1,:), 'r', 'DisplayName', ['DeePC y1, lg=', num2str(l_g)]);
stairs(t, y(2,:), 'b', 'DisplayName', 'DeePC y2');
ylabel('Output [-]');
legend;
grid on; grid minor;
axis([0 t(end) -1 15]);


ax3 = subplot(212);
hold on;
stairs(t, u(1,:), 'r');
stairs(t, u(2,:), 'b');
%stairs(t, u_ss(2)*ones(1,length(t)), '--k');
axis([0 t(end) -20 20])
ylabel('Input [-]');
xlabel('Time [s]');
grid on; grid minor;

%%

CtrlCost = 0;
for k = Tini+1:simLen
    CtrlCost = CtrlCost + (y(:,k)-r(:,k))'*Q*(y(:,k)-r(:,k))+u(:,k)'*R*u(:,k);
end
disp(['DeePC control cost =', num2str(CtrlCost)]);

