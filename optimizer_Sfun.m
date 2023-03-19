function [sys,x0,str,ts] = optimizer_Sfun(t,z,y,flag)



switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes();

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,z,y);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,z,y);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes()


sizes = simsizes;
sizes.NumContStates  = 22;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);


x0  = rand(1,22)*100;
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,z,y,Parameters)


load("Parameters/optimizer_parameters.mat")
load("Parameters/growth_parameters.mat")



% Optimizer
xd  = z(1);
yd  = z(2);
o1  = z(3);
o2  = z(4);
r   = z(5);
a   = z(6);
xm  = z(7);  % gRNA
xp  = z(8);  % gRNA
ym  = z(9);  % gRNA
yp  = z(10); % gRNA
xmb = z(11); % STAR
xpb = z(12); % STAR
ymb = z(13); % STAR
ypb = z(14); % STAR
qpp = z(15); 
qmm = z(16);
qpm = z(17);
qmp = z(18);
u1  = z(19);
u2  = z(20);
x   = z(21);
v   = z(22);

% Parameters of the optimizer
delta_c  = par_opt.delta_c;
delta_q  = par_opt.delta_q;

lambda_d = par_opt.lambda_d;
lambda_o = par_opt.lambda_o;
lambda_u = par_opt.lambda_u;
lambda_x = par_opt.lambda_x;
lambda_v = par_opt.lambda_v;

alpha_x  = par_opt.alpha_x;
alpha_y  = par_opt.alpha_y;
alpha_xb = par_opt.alpha_xb;
alpha_yb = par_opt.alpha_yb;
alpha_q  = par_opt.alpha_q;

beta_x   = par_opt.beta_x;
beta_xb  = par_opt.beta_xb;
beta_v   = par_opt.beta_v;
beta_xd  = par_opt.beta_xd;
beta_yd  = par_opt.beta_yd;
beta_u   = par_opt.beta_u;
beta_o   = par_opt.beta_o;

k_xb     = par_opt.k_xb;
k_yb     = par_opt.k_yb;
k_xt     = par_opt.k_xt;
k_yt     = par_opt.k_yt;
k_q      = par_opt.k_q;

K_o      = par_opt.K_o;
K_x      = par_opt.K_x;
K_u1     = par_opt.K_u1;
K_u2     = par_opt.K_u2;
K_y      = par_opt.K_y;
K_xb     = par_opt.K_xb;
K_yb     = par_opt.K_yb;
K_c      = par_opt.K_c;

nu       = par_opt.nu;

% Parameters that are shared: growth rate
mu       = par_growth.mu;



% Delay    
sys(1)  = beta_xd*x/(x+K_x) - (mu + lambda_d)*xd;
sys(2)  = beta_yd*y/(y+K_y) - (mu + lambda_d)*yd;

% Repressilator
sys(3)  = beta_o/(1+(r/K_o)^4) - (mu+lambda_o)*o1;
sys(4)  = beta_o/(1+(o1/K_o)^2) - (mu+lambda_o)*o2;
sys(5)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*r;
sys(6)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*a;

% gRNAs
sys(7)  = alpha_x*K_c^4/(r^4+K_c^4)*xd/(xd+K_xb) + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xm;
sys(8)  = alpha_x*K_c^4/(r^4+K_c^4)*x/(x+K_xb)   + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xp;
sys(9)  = alpha_y*K_c^4/(r^4+K_c^4)*yd/(yd+K_yb) + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ym;
sys(10) = alpha_y*K_c^4/(r^4+K_c^4)*y/(y+K_yb)   + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*yp;

% STARs
sys(11)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xmb;
sys(12)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xpb;
sys(13)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ymb;
sys(14)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*ypb;

% AND gates
sys(15) = alpha_q*xpb/(xpb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qpp;
sys(16) = alpha_q*xmb/(xmb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qmm;
sys(17) = alpha_q*xpb/(xpb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qpm;
sys(18) = alpha_q*xmb/(xmb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qmp;

% OR gates
sys(19) = beta_u*(qpp/(k_q + qpp)+qmm/(k_q + qmm)) - (mu + lambda_u)*u1;
sys(20) = beta_u*(qpm/(k_q + qpm)+qmp/(k_q + qmp)) - (mu + lambda_u)*u2;

% Regulator
sys(21) = beta_x*x/(x+K_x) + beta_xb*u1^2/(u1^2+K_u1^2) - (mu + lambda_x*(1+nu*v))*x;

% Protease
sys(22) = beta_v*u2/(u2+K_u2) - (mu + lambda_v)*v;



% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,z,y)


sys = z(21);


% end mdlOutputs