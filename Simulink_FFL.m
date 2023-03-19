clear all
close all
clc

% simulation horizon
T_end = 300;

% temporal pattern of w
delta = 1e-6;
t = [0 T_end/3-delta T_end/3+delta 2*T_end/3-delta 2*T_end/3+delta T_end];
w = [0 0 10 10 5 5];
w_t = [t' w'];

% load parameters
load('Parameters/growth_parameters.mat')
load('Parameters/process_parameters.mat')

% optimum of x
x_opt = par_proc.kappa_x*sqrt(par_proc.kappa_z*(w+par_proc.kappa_w)/par_proc.kappa_w*(par_growth.mu+par_proc.lambda_z)/par_proc.beta_z);

% simulate closed loop system
sim_data = sim('ClosedLoop_FFL.slx');
sim_t    = sim_data.simout.Time;
sim_x    = sim_data.simout.Data;

% plot the results
figure()
hold on
plot(t,x_opt,'b','LineWidth',5)      % plot optimum
plot(sim_t,sim_x,'r','LineWidth',5)  % plot closed loop behavior

% styling figure
ylim([0 100]), grid on
xlabel('time (h)'), ylabel('x (nM)')
legend('optimum','closed loop')
