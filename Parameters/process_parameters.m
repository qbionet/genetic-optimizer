% Parameters of the process to be regulated
% Only uncomment the one that is needed, comment the rest
clear all
close all
clc


% Simple reporter
% par_proc.lambda_y = 1;      % degradation rate constant of the reporter (1/h)
% par_proc.x_opt = 30;        % location of optimum (nM)
% par_proc.F0 = 100;          % height of objective function (nM/h)
% par_proc.sigma = 10;        % width of objective function (nM)
% par_proc.dim   = 1;         % dimension of the process (-)


%% FFL 
par_proc.lambda_y = 1;        % degradation rate constant of y (1/h)
par_proc.lambda_z = 2;        % degradation rate constant of z (1/h)
par_proc.beta_y   = 50;       % production rate constant of y (nM/h)
par_proc.beta_z   = 1200;     % production rate constant of z (nM/h)
par_proc.kappa_x  = 70;       % dissociation rate constant of x (nM)
par_proc.kappa_z  = 200;      % dissociation rate constant of z (nM)
par_proc.kappa_w  = 10;       % dissociation rate constant of w (nM)
par_proc.dim   = 2;           % dimension of the process (-)


%% save parameters
save("process_parameters.mat")

