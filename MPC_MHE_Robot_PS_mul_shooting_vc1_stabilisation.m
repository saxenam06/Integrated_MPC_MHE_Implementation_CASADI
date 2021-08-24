% First CASADI test for coupled MPC and MHE for Mobile Robots
clear all
close all
clc

addpath('C:\CASADI')
import casadi.*

T = 0.1; % Time Step [s]
N = 30; % prediction horizon
N_MHE = 30; % Estimation horizon
rob_diam = 0.3;

%%%% System Model
x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');  % States
states = [x;y;theta]; n_states = length(states);    
v = SX.sym('v'); omega = SX.sym('omega');   % Controls
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % r.h.s of MOTION MODEL
f = Function('f',{states,controls},{rhs}); % System Model

%----------------------MPC Problem  setup starts%%%%%%%%%%%%%%%%%%%%
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);% parameters (which include initial state of the robot and the reference state)
X = SX.sym('X',n_states,(N+1));% A vector that represents the states over the prediction horizon.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.5; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.05; R(2,2) = 0.05; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % Equality constraints defined by Discrete motion model
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);
args = struct;

args.lbg(1:3*(N+1)) = 0;
args.ubg(1:3*(N+1)) = 0;

args.lbx(1:3:3*(N+1),1) = -4; %state x lower bound
args.ubx(1:3:3*(N+1),1) = 4; %state x upper bound
args.lbx(2:3:3*(N+1),1) = -4; %state y lower bound
args.ubx(2:3:3*(N+1),1) = 4; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) = inf; %state theta upper bound

v_max = 0.6; v_min = -v_max;
omega_max = pi/4; omega_min = -omega_max;

args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %v lower bound
args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %v upper bound
args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_min; %omega lower bound
args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_max; %omega upper bound
%----------------------MPC Problem setup done%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%MHE Problem setup starts%%%%%%%%%%%%%%%%%%%
% MEASUREMENT MODEL
r = SX.sym('r'); alpha = SX.sym('alpha'); % range and bearing
measurement_rhs = [sqrt(x^2+y^2); atan(y/x)];
h = Function('h',{states},{measurement_rhs}); % MEASUREMENT MODEL

U_est = SX.sym('U_est',n_controls,N_MHE);    %(controls)
X_est = SX.sym('X_est',n_states,(N_MHE+1));  %(states) 
P_meas = SX.sym('P_meas', 2 , (N_MHE+1)+N_MHE); % parameters (include range & bearing sensor measurements and calculated controls)

con_cov = diag([0.001 deg2rad(1.5)]).^2;
meas_cov = diag([0.001 deg2rad(1.5)]).^2;

V = inv(sqrt(meas_cov)); % weighing matrices (output)  y_tilde - y
W = inv(sqrt(con_cov)); % weighing matrices (input)   u_tilde - u

obj = 0; % Objective function
g = [];  % constraints vector
for k = 1:N_MHE+1
    st = X_est(:,k);
    h_x = h(st);
    y_tilde = P_meas(:,k);
    obj = obj+ (y_tilde-h_x)' * V * (y_tilde-h_x); % calculate obj
end

for k = 1:N_MHE
    con = U_est(:,k);
    u_tilde = P_meas(:,N_MHE+1+k);
    obj = obj+ (u_tilde-con)' * W * (u_tilde-con); % calculate obj
end

% multiple shooting constraints
for k = 1:N_MHE
    st = X_est(:,k);  con = U_est(:,k);
    st_next = X_est(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % Equality constraints defined by Discrete motion model
end


% make the decision variable one column  vector
OPT_variables_est = [reshape(X_est,3*(N_MHE+1),1);reshape(U_est,2*N_MHE,1)];

nlp_mhe = struct('f', obj, 'x', OPT_variables_est, 'g', g, 'p', P_meas);

solver_mhe = nlpsol('solver', 'ipopt', nlp_mhe,opts);


args1 = struct;

args1.lbg(1:3*(N_MHE)) = 0;
args1.ubg(1:3*(N_MHE)) = 0;

args1.lbx(1:3:3*(N_MHE+1),1) = -4; %state x lower bound
args1.ubx(1:3:3*(N_MHE+1),1) = 4; %state x upper bound
args1.lbx(2:3:3*(N_MHE+1),1) = -4; %state y lower bound
args1.ubx(2:3:3*(N_MHE+1),1) = 4; %state y upper bound
args1.lbx(3:3:3*(N_MHE+1),1) = -pi/2; %state theta lower bound
args1.ubx(3:3:3*(N_MHE+1),1) = pi/2; %state theta upper bound

args1.lbx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_min; %v lower bound
args1.ubx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_max; %v upper bound
args1.lbx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_min; %omega lower bound
args1.ubx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_max; %omega upper bound
%----------------------------------------------

%%%%%%%%%%%%%%%%%%%%%MHE Problem setup done%%%%%%%%%%%%%%%%%%%


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
t(1) = t0;
sim_tim = 40; % total sampling times

x0 = [0 ; 0 ; 0];    % initial condition.
xg = [0 ; 0 ; 0];    % initial condition.
xs = [3.5 ; 3.5 ; 0.0]; % Reference posture.

%%% Simulating first measurement at initial condition, No noise
h_value0=full(h(x0));  
y_measurements = [h_value0(1,1) h_value0(2,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X0_est = zeros(N_MHE+1,3); % used for initialization of the MHE estimated states decision variables as per the measurements later (|)
U0_est = zeros(N_MHE,2);   % used for initialization of the MHE estimated two control inputs as per the planned control inputs (|)

X0 = repmat(x0,1,N+1)'; % used for initialization of the MPC states decision variables (|)
u0 = zeros(N,2);        % used initialisation of the MPC control inputs decision variables (|)

% Start MPC
mpciter = 0;
mheiter = 0;

xx_G(:,1) = xg; % contains the history of ground truth for MHE(-)
xx_m(:,1) = x0(1:2,1);  % contains the history of x & y directly calculated from sensor measurements (-)

xx(:,1) = x0; % xx contains the history of exhibited states for MPC(-)
xx1 = [];   % for storing predicted trajectories from MPC
u_cl=[];    % for storing planned measured control action from MPC (|)

tic;
while(norm((x0-xs),2) > 0.005 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution (|)
    xx1(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY (|)
    u_cl= [u_cl ; u(1,:)]; %(|)
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0,r ,alpha,rn ,alphan,u0,x0,xg ] = shift_meas_control_noise(T, t0, x0, u,f,h,xg);
    % Storing Ground Truth
    xx_G(:,mpciter+2) = xg;
    % Collecting Measurements until N_MHE prediction horizon as per ground
    % truth
    if mpciter+2 <= N
        %%%if not enough measurement samples yet equal to MHE Prediction
        %%% horizon assume ground truth known
    y_measurements = [y_measurements; r , alpha];
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % Initialise the states decision variables for MPC with ground truth until MHE est not available
    X0 = [X0(2:end,:);X0(end,:)];  
    xx(:,mpciter+2) = x0; %% Assuming ground truth is known until NMHE steps
    xx_m(:,mpciter+2) = x0(1:2,1); % %store of x & y directly calculated from sensor measurements (-)
    else
        %%if enough samples  equal to MHE Prediction
        %%% horizon then perform state estimation using MHE
        y_measurements = [y_measurements; rn , alphan];
        xx_m(:,mpciter+2) = [rn*cos(alphan),rn*sin(alphan)]'; % storing measurements x & y directly calculated from sensor measurements (-)
        if mheiter==0
            U0_est = u_cl(1:N_MHE,:); % initialize the control actions by the earlier calculated control law
            X0_est(:,1:3) = [y_measurements(1:N_MHE+1,1).*cos(y_measurements(1:N_MHE+1,2)),...
            y_measurements(1:N_MHE+1,1).*sin(y_measurements(1:N_MHE+1,2)),[xx(3,:) xx(3,end)]']; % initialize the states from the measured range and bearing
            mheiter = mheiter + 1
        else
            X0_est = [X_sol(2:end,:);X_sol(end,:)];
            U0_est = [U_sol(2:end,:);U_sol(end,:)];
            mheiter = mheiter + 1
        end
        % Get the measurements window and put it as parameters in p
        args1.p   = [y_measurements(mheiter:mheiter+N_MHE,:)',u_cl(mheiter:mheiter+N_MHE-1,:)'];
        % initial value of the optimization variables
        args1.x0  = [reshape(X0_est',3*(N_MHE+1),1);reshape(U0_est',2*N_MHE,1)];
        sol_mhe = solver_mhe('x0', args1.x0, 'lbx', args1.lbx, 'ubx', args1.ubx,...
        'lbg', args1.lbg, 'ubg', args1.ubg,'p',args1.p);
        U_sol = reshape(full(sol_mhe.x(3*(N_MHE+1)+1:end))',2,N_MHE)'; % get controls only from the solution
        X_sol = reshape(full(sol_mhe.x(1:3*(N_MHE+1)))',3,N_MHE+1)'; % get solution TRAJECTORY
        x0=X_sol(N_MHE+1,:)';
        xx(:,mpciter+2) =x0;
    end
    mpciter;
    mpciter = mpciter + 1
end
toc

ss_error = norm((x0-xs),2)
Draw_MPC_MHE_point_stabilization_v1(t,xx_G,xx_m,xx,xx1,u_cl,xs,N,rob_diam)

figure(1)
set(gca,'Fontsize',14)
subplot(311)
plot(t,xx_G(1,1:end-1),'-.b','linewidth',1.5); hold on
ylabel('x (m)')
grid on
set(gca,'Fontsize',14)
subplot(312)
plot(t,xx_G(2,1:end-1),'-.b','linewidth',1.5);hold on
ylabel('y (m)')
grid on
set(gca,'Fontsize',14)
subplot(313)
plot(t,xx_G(3,1:end-1),'-.b','linewidth',1.5); hold on
xlabel('time (seconds)')
ylabel('\theta (rad)')
grid on

% Plot the cartesian coordinates from the measurements used
figure(1)
set(gca,'Fontsize',14)
subplot(311)
plot(t,y_measurements(1:end-1,1).*cos(y_measurements(1:end-1,2)),'.m','linewidth',1.5); hold on
grid on
set(gca,'Fontsize',14)
subplot(312)
plot(t,y_measurements(1:end-1,1).*sin(y_measurements(1:end-1,2)),'.m','linewidth',1.5); hold on
grid on


figure(1)
set(gca,'Fontsize',14)
subplot(311)
plot(t,xx(1,1:end-1),'--g','linewidth',1.5); hold on
ylabel('x (m)')
legend('Ground Truth','Measurement','Exhibited Trajectory')
grid on
set(gca,'Fontsize',14)
subplot(312)
plot(t,xx(2,1:end-1),'--g','linewidth',1.5);hold on
ylabel('y (m)')
legend('Ground Truth','Measurement','Exhibited Trajectory')
grid on
set(gca,'Fontsize',14)
subplot(313)
plot(t,xx(3,1:end-1),'--g','linewidth',1.5); hold on
legend('Ground Truth','Exhibited Trajectory')
xlabel('time (seconds)')
ylabel('\theta (rad)')
grid on
print -dpng 'states_exhibit_meas'



% plot the ground truth mesurements VS the noisy measurements
figure(2)
% set(y1, 'PaperUnits','centimeters','PaperPosition',[0 0 22 18] )
set(gca,'Fontsize',14)
subplot(211)
plot(t,sqrt(xx(1,1:end-1).^2+xx(2,1:end-1).^2),'b','linewidth',1.5); hold on
plot(t,y_measurements(1:end-1,1),'r','linewidth',1.5); hold on
ylabel('Range: [ r (m) ]')
grid on
legend('Ground Truth','Measurement')
set(gca,'Fontsize',14)
subplot(212)
plot(t,atan(xx(2,1:end-1)./xx(1,1:end-1)),'b','linewidth',1.5); hold on
plot(t,y_measurements(1:end-1,2),'r','linewidth',1.5);hold on
ylabel('Bearing: [ \alpha (rad) ]')
print -dpng 'Sensor_noise'
grid on

