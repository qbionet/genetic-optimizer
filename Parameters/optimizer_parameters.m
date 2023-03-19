clear all
close all
clc


% Parameters of the optimizer (Supplementary Table 1)
par_opt.delta_c  = 20;      % gRNA/STAR degradation rate constant in comparator (1/h)
par_opt.delta_q  = 6;       % STAR degradation rate constant in AND gates (1/h)

par_opt.lambda_x = 0.2;     % degradation rate constant of the regulator (1/h)
par_opt.lambda_v = 4;       % degradation rate constant of the protease (1/h)
par_opt.lambda_d = 4;       % degradation rate constant of the delayed species (1/h)
par_opt.lambda_u = 6;       % degradation rate constant of the control signals (1/h)
par_opt.lambda_o = 6;       % degradation rate constant of the repressors in the oscillator (1/h)

par_opt.alpha_x  = 21e3;    % RNA production rate constant during phase 1 in the comparator (nM/h)
par_opt.alpha_y  = 2.1e3;   % RNA production rate constant during phase 1 in the comparator (nM/h)
par_opt.alpha_xb = 21e3;    % RNA production rate constant during phase 2 in the comparator (nM/h)
par_opt.alpha_yb = 10.5e3;  % RNA production rate constant during phase 2 in the comparator (nM/h)
par_opt.alpha_q  = 0.7e3;   % RNA production rate constant in the AND gates (nM/h)

par_opt.beta_x   = 6.1e3;   % protein production rate constant of the regulator due to self-activation (nM/h)
par_opt.beta_xb  = 0.12e3;  % protein production rate constant of the regulator due to regulation (nM/h)
par_opt.beta_v   = 2e3;     % protein production rate constant of the protease (nM/h)
par_opt.beta_xd  = 25e3;    % protein production rate constant of the delayed signal of the regulator (nM/h)
par_opt.beta_yd  = 25e3;    % protein production rate constant of the delayed signal of the reporter (nM/h)
par_opt.beta_u   = 3.5e3;   % protein production rate constant of control signals (nM/h)
par_opt.beta_o   = 20e3;    % protein production rate constant of the repressors in the oscillator (nM/h)

par_opt.k_xb     = 100;     % dissociation rate constant of gRNA with terminator hairpin (nM)
par_opt.k_yb     = 10;      % dissociation rate constant of gRNA with terminator hairpin (nM)
par_opt.k_xt     = 1e3;     % dissociation rate constant of STAR with terminator hairpin (nM)
par_opt.k_yt     = 1e3;     % dissociation rate constant of STAR with terminator hairpin (nM)
par_opt.k_q      = 500;     % dissociation rate constant of STAR with terminator hairpin (nM)

par_opt.K_o      = 1;       % dissociation rate constant of repressors in repressilator (nM)
par_opt.K_c      = 10;      % dissociation rate constant of repressors in delay modules (nM)
par_opt.K_x      = 5e3;     % dissociation rate constant of regulator self-activation (nM)
par_opt.K_y      = 5e3;     % dissociation rate constant of reporter with cognate promoter (nM)
par_opt.K_u1     = 40;      % dissociation rate constant of u1 with cognate promoter (nM)
par_opt.K_u2     = 200;     % dissociation rate constant of u2 with cognate promoter (nM)
par_opt.K_xb     = 1e3;     % dissociation rate constant of x and x_d in comparator (nM)
par_opt.K_yb     = 100;     % dissociation rate constant of y and y_d in comparator (nM)

par_opt.nu       = 0.05;    % parameter affecting the controlled degradation of the regulator (1/nM)

save('optimizer_parameters.mat')


