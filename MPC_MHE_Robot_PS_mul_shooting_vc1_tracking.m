% First CASADI test for coupled MPC and MHE for Mobile Robots
clear all
close all
clc

addpath('C:\CASADI')
import casadi.*

T = 0.1; %[s]
N = 50; % prediction horizon
N_MHE = 50; % Estimation horizon
rob_diam = 0.3;
xs = [12 ; 1 ; 0.0]; % Reference posture.

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s MOTION MODEL

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + N*(n_states+n_controls));
X = SX.sym('X',n_states,(N+1));% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.75; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.05; R(2,2) = 0.05; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(5*k-1:5*k+1))'*Q*(st-P(5*k-1:5*k+1)) + ...
              (con-P(5*k+2:5*k+3))'*R*(con-P(5*k+2:5*k+3)) ; % calculate obj    st_next = X(:,k+1);
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
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

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

args.lbg(1:3*(N+1)) = 0;
args.ubg(1:3*(N+1)) = 0;

args.lbx(1:3:3*(N+1),1) = -20; %state x lower bound
args.ubx(1:3:3*(N+1),1) = 20; %state x upper bound
args.lbx(2:3:3*(N+1),1) = -2; %state y lower bound
args.ubx(2:3:3*(N+1),1) = 2; %state y upper bound
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
r = SX.sym('r'); alpha = SX.sym('alpha'); % range and bearing
measurement_rhs = [sqrt(x^2+y^2); atan(y/x)];
h = Function('h',{states},{measurement_rhs}); % MEASUREMENT MODEL

U_est = SX.sym('U_est',n_controls,N_MHE);    %(controls)
X_est = SX.sym('X_est',n_states,(N_MHE+1));  %(states) [remember multiple shooting]
P_meas = SX.sym('P_meas', 2 , (N_MHE+1)+N_MHE); % parameters (include r and alpha measurements as well as controls measurements)

con_cov = diag([0.001 deg2rad(0.2)]).^2;
% add noise to the range and bearing sensors
meas_cov = diag([0.001 deg2rad(0.2)]).^2;
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
    g = [g;st_next-st_next_euler]; % compute constraints
end


% make the decision variable one column  vector
OPT_variables_est = [reshape(X_est,3*(N_MHE+1),1);reshape(U_est,2*N_MHE,1)];

nlp_mhe = struct('f', obj, 'x', OPT_variables_est, 'g', g, 'p', P_meas);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver_mhe = nlpsol('solver', 'ipopt', nlp_mhe,opts);


args1 = struct;

args1.lbg(1:3*(N_MHE)) = 0;
args1.ubg(1:3*(N_MHE)) = 0;

args1.lbx(1:3:3*(N_MHE+1),1) = -20; %state x lower bound
args1.ubx(1:3:3*(N_MHE+1),1) = 20; %state x upper bound
args1.lbx(2:3:3*(N_MHE+1),1) = -2; %state y lower bound
args1.ubx(2:3:3*(N_MHE+1),1) = 2; %state y upper bound
args1.lbx(3:3:3*(N_MHE+1),1) = -pi/2; %state theta lower bound
args1.ubx(3:3:3*(N_MHE+1),1) = pi/2; %state theta upper bound

args1.lbx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_min; %v lower bound
args1.ubx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_max; %v upper bound
args1.lbx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_min; %omega lower bound
args1.ubx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_max; %omega upper bound
%----------------------------------------------

% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
t(1) = t0;
sim_tim = 40; % total sampling times

x0 = [0.5 ; 0.5 ; 0];    % initial condition.
xg = [0.5 ; 0.5 ; 0];    % initial condition.
xx_G(:,1) = xg; % xx contains the history of ground truth (-)
xx_m(:,1) = x0(1:2,1);
h_value0=full(h(x0));
xx(:,1) = x0; % xx contains the history of estimated states (-)
xx1 = [];   % for storing predicted trajectories from MPC

u_cl=[];    % for storing planned measured control action from MPC (|)
y_measurements = [h_value0(1,1) h_value0(2,1)];

X0_est = zeros(N_MHE+1,3); % initialization of the states decision variables as per the measurements later (|)
U0_est = zeros(N_MHE,2);   % two control inputs for each robot as per the planned control inputs (|)

X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables (|)
u0 = zeros(N,2);        % initialisation of the control inputs decision variables (|) 

% Start MPC
mpciter = 0;
mheiter = 0;

% r=[];
% alpha=[];
% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
% tic
while(norm((x0-xs),2) > 0.5 && mpciter < sim_tim / T)
    current_time = mpciter*T;  %new - get the current time
    % args.p   = [x0;xs]; % set the values of the parameters vector
    %----------------------------------------------------------------------
    args.p(1:3) = x0; % initial condition of the robot posture
    for k = 1:N %new - set the reference to track
        t_predict = current_time + (k-1)*T; % predicted time instant
        x_ref = 0.5*t_predict; y_ref = 1; theta_ref = 0;
        u_ref = 0.5; omega_ref = 0;
        if x_ref >= 12 % the trajectory end is reached
            x_ref = 12; y_ref = 1; theta_ref = 0;
            u_ref = 0; omega_ref = 0;
        end
        args.p(5*k-1:5*k+1) = [x_ref, y_ref, theta_ref];
        args.p(5*k+2:5*k+3) = [u_ref, omega_ref];
    end
    % initial value of the optimization variables
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution (|)
    xx1(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY (|)
    u_cl= [u_cl ; u(1,:)]; %(|)
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0,r ,alpha,rn ,alphan,u0,x0,xg ] = shift_meas_control_noise_track(T, t0, x0, u,f,h,xg);
    xx_G(:,mpciter+2) = xg;
    if mpciter+2 <= N
    y_measurements = [y_measurements; r , alpha];
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    X0 = [X0(2:end,:);X0(end,:)];  
    xx(:,mpciter+2) = x0; %%Assuming ground truth is known until NMHE steps
    xx_m(:,mpciter+2) = x0(1:2,1); % storing measurements
    else
    y_measurements = [y_measurements; rn , alphan];
    xx_m(:,mpciter+2) = [rn*cos(alphan),rn*sin(alphan)]'; % storing measurements
        if mheiter==0
            U0_est = u_cl(1:N_MHE,:); % initialize the control actions by the measured
            X0_est(:,1:3) = [y_measurements(1:N_MHE+1,1).*cos(y_measurements(1:N_MHE+1,2)),...
                y_measurements(1:N_MHE+1,1).*sin(y_measurements(1:N_MHE+1,2)),[xx(3,:) xx(3,end)]']; % initialize the states from the measured range and bearing
            mheiter = mheiter + 1
            
        else
            X0_est = [X_sol(2:end,:);X_sol(end,:)];
            U0_est = [U_sol(2:end,:);U_sol(end,:)];
            mheiter = mheiter + 1
        end
%         tic
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
        % Shift trajectory to initialize the next step
%         toc
    
    mpciter;
    mpciter = mpciter + 1

end
% toc

Draw_MPC_tracking_v1(t,xx_G,xx_m,xx,xx1,u_cl,N,rob_diam)


%-----------------------------------------
%-----------------------------------------
%-----------------------------------------
%    Start MHE implementation from here
%-----------------------------------------
%-----------------------------------------
%-----------------------------------------
% plot the ground truth
figure(1)
subplot(311)
plot(t,xx_G(1,1:end-1),'-.b','linewidth',1.5); hold on
ylabel('x (m)')
grid on
subplot(312)
plot(t,xx_G(2,1:end-1),'-.b','linewidth',1.5);hold on
ylabel('y (m)')
grid on
subplot(313)
plot(t,xx_G(3,1:end-1),'-.b','linewidth',1.5); hold on
xlabel('time (seconds)')
ylabel('\theta (rad)')
grid on

% Plot the cartesian coordinates from the measurements used
figure(1)
subplot(311)
plot(t,y_measurements(1:end-1,1).*cos(y_measurements(1:end-1,2)),'r','linewidth',1.5); hold on
grid on
legend('Ground Truth','Measurement')
subplot(312)
plot(t,y_measurements(1:end-1,1).*sin(y_measurements(1:end-1,2)),'r','linewidth',1.5); hold on
grid on

% plot the ground truth mesurements VS the noisy measurements
figure(2)
subplot(211)
plot(t,sqrt(xx(1,1:end-1).^2+xx(2,1:end-1).^2),'b','linewidth',1.5); hold on
plot(t,y_measurements(1:end-1,1),'r','linewidth',1.5); hold on
ylabel('Range: [ r (m) ]')
grid on
legend('Ground Truth','Measurement')
subplot(212)
plot(t,atan(xx(2,1:end-1)./xx(1,1:end-1)),'b','linewidth',1.5); hold on
plot(t,y_measurements(1:end-1,2),'r','linewidth',1.5);hold on
ylabel('Bearing: [ \alpha (rad) ]')
grid on





figure(1)
subplot(311)
plot(t,xx(1,1:end-1),'--g','linewidth',1.5); hold on
ylabel('x (m)')
grid on
subplot(312)
plot(t,xx(2,1:end-1),'--g','linewidth',1.5);hold on
ylabel('y (m)')
grid on
subplot(313)
plot(t,xx(3,1:end-1),'--g','linewidth',1.5); hold on
xlabel('time (seconds)')
ylabel('\theta (rad)')
grid on

