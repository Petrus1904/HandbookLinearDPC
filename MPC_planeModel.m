%% Intro
% MPC simulation of a Plane model. Considers classic MPC control law
% without a stabilizing terminal cost/constraint. The system dynamics are obtained from
% a data set and identified using N4SID. The state is observed using a
% Luenberger Observer, where the observer gain is the result of the LQR
% dual problem. Data acquisition and simulation is subject to sensor noise

clear all;

%Set this to any constant to fix the noise seed every run.
rng('shuffle');

%% Define controller parameters and obtain Hankel matrices
Ts = 0.1;
N = 20;
Tini = 20;
T = 2500;

wvar = 0.25; %noise variance

[sys, constr, Hankel] = GetDataPlaneModel(T, N, Tini, wvar, Ts);

A = sys.A;
B = sys.B;
C = sys.C;

n = size(A,1);
nu = size(B,2);
ny = size(C,1);

u_ss = inv([A-eye(n), B; C, zeros(ny,nu)])*[zeros(n,1); [10;0]];
u_ss = u_ss(end-nu+1:end);

%n4sid uses a random call, thus ruining the RNG of the monte carlo runs.
%Hence I store it here
RS = rng;
Iddat = iddata(Hankel.Y', Hankel.U', Ts);
SYSid = n4sid(Iddat,4); %we assume perfect prediction, despite n4sid indicating n=3 is a better fit
%and return to MY random sequence here.
rng(RS);

Cid = SYSid.C;
Bid = SYSid.B;
Aid = SYSid.A;
nid = size(Aid,1);
Kid = dlqr(Aid', Cid', eye(nid), zeros(ny))';
Phi = GetMarkovMatrix(Cid, Aid, eye(nid), (1:N)');
Gamma = GetMarkovMatrix(Cid, Aid, Bid, tril(ones(N))*tril(ones(N))-1);

%% Build controller
Q = 10*eye(ny);
R = 0.01*eye(nu);

Psi = kron(eye(N), R);
Omega = kron(eye(N), Q);

cpu=[];
u = sdpvar(nu*N,1);
y = sdpvar(ny*N,1);
ref = sdpvar(ny*N,1);
xh = sdpvar(nid,1);

% Define objective function
objective = (y-ref)'*Omega*(y-ref)+(u)'*Psi*(u); % base MPC cost

% Build constraints
constraints = [y == Phi*xh + Gamma*u];

for k = 1:N   % and the bounds on the inputs and outputs
    constraints = [constraints, constr.umin<=u(nu*(k-1)+1:nu*k)<=constr.umax, constr.ymin<=y(ny*(k-1)+1:ny*k)<=constr.ymax];
end

Parameters = {xh, ref};
Outputs = {u, y};
options = sdpsettings('solver', 'osqp', 'verbose', 0, 'debug', 0); % 'osqp.eps_abs', 1e-8, 'osqp.eps_rel', 1e-8);
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


rng(1904);
y = zeros(ny, simLen);
u = zeros(nu, simLen);
x = zeros(n, simLen+1);
x(:,1) = x0;
xh = zeros(nid, simLen+1);

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

        [Sol, err] = controller({xh(:,k), Rk});
        if(err ~= 0)
            stop = 1;
        end
        Uk = Sol{1}; 
        Yk = Sol{2}; %Yk can be used for evaluation, is not stored
        cpu = [cpu toc];
       
        u(:,k) = Uk(1:nu);
        xh(:,k+1) = Aid*xh(:,k) + Bid*u(:,k)+Kid*(y(:,k) - Cid*xh(:,k));
        
    else
        %even though model based controller do not need an initial window,
        %we did it for consistency with the DPC laws
        
        u(:,k) = 1*randn(nu, 1);
        xh(:,k+1) = Aid*xh(:,k) + Bid*u(:,k)+Kid*(y(:,k) - Cid*xh(:,k));
    end

    % Update system
    x(:,k+1) = A*x(:,k) + B*u(:,k);
end

%% Display simulation results

% figure();
% ax1 = subplot(211);
% stairs(t, r(1,1:simLen), 'k', 'DisplayName', 'Reference 1');
% hold on;
% stairs(t, r(2,1:simLen), 'k', 'DisplayName', 'Reference 2');
% stairs(t, y(1,:), 'r', 'DisplayName', 'MPC y1');
% stairs(t, y(2,:), 'b', 'DisplayName', 'MPC y2');
% stairs(t, Cid(1,:)*xh(:,1:end-1), 'g--', 'DisplayName', 'Obs y1');
% stairs(t, Cid(2,:)*xh(:,1:end-1), 'm--', 'DisplayName', 'Obs y2');
% ylabel('Output [-]');
% legend;
% grid on; grid minor;
% axis([0 t(end) -1 15]);
% 
% 
% ax3 = subplot(212);
% hold on;
% stairs(t, u(1,:), 'r');
% stairs(t, u(2,:), 'b');
% %stairs(t, u_ss(2)*ones(1,length(t)), '--k');
% axis([0 t(end) -20 20])
% ylabel('Input [-]');
% xlabel('Time [s]');
% grid on; grid minor;

%%

CtrlCost = 0;
for k = Tini+1:simLen
    CtrlCost = CtrlCost + (y(:,k)-r(:,k))'*Q*(y(:,k)-r(:,k))+u(:,k)'*R*u(:,k);
end
disp(['MPC control cost =', num2str(CtrlCost)]);
