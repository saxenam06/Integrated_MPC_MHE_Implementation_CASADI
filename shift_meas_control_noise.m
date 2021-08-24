function [t0,r ,alpha,rn ,alphan,u0,x0,xg ] = shift_meas_control_noise(T, t0, x0, u,f,h,xg)
% add noise to the control actions 
con_cov = diag([0.001 deg2rad(1.5)]).^2;
% add noise to the range and bearing sensors
meas_cov = diag([0.001 deg2rad(1.5)]).^2;


con = u(1,:)'+ sqrt(con_cov)*randn(2,1); 


st = x0;
f_value = f(st,con);   
st = st+ (T*f_value);
x0 = full(st);

%%%%%Simulating Ground Truth
st = xg;
f_value = f(st,con);   
st = st+ (T*f_value);
xg = full(st);

% %%%%%Simulating Sensor measurement by adding noise on the Ground truth

h_value=full(h(xg));
rn = h_value(1,1)+sqrt(meas_cov(1,1))*randn(1);
alphan = h_value(2,1)+sqrt(meas_cov(2,2))*randn(1);

r = h_value(1,1);
alpha = h_value(2,1);

t0 = t0 + T;
u0 = [u(2:size(u,1),:);u(size(u,1),:)]; % shift the control action 
end