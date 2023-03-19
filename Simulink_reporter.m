clear all
close all
clc

% simulation horizon
T_end = 300;

% temporal pattern of x_opt
delta   = 1e-6;
t       = [0 T_end/3-delta T_end/3+delta 2*T_end/3-delta 2*T_end/3+delta T_end];
x_opt   = [20 20 80 80 50 50];
x_opt_t = [t' x_opt'];

% load parameters
load('Parameters/growth_parameters.mat')
load('Parameters/process_parameters.mat')


% simulate closed loop system
sim_data = sim('ClosedLoop_reporter.slx');
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
