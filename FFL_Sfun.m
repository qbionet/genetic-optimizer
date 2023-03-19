function [sys,x0,str,ts] = process_Sfun(t,state,input,flag)



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
    sys=mdlDerivatives(t,state,input);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,state,input);

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
sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);


x0  = rand(1,2)*100;
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,state,input,Parameters)

load("Parameters/process_parameters.mat")
load("Parameters/growth_parameters.mat")

% Parameters of the process


% Parameters that are shared: growth rate
mu       = par_growth.mu;       % growth rate (1/h)



% Parameters of the process (these need to be defined in process_parameters.m)
lambda_y  = par_proc.lambda_y;

beta_z    = par_proc.beta_z;
beta_y    = par_proc.beta_y;
lambda_z  = par_proc.lambda_z;
kappa_x   = par_proc.kappa_x;
kappa_z   = par_proc.kappa_z;
kappa_w   = par_proc.kappa_w;

z = state(1);
y = state(2);

x = input(1);
w = input(2);

sys(1) = beta_z*(x/kappa_x)/(1 + w/kappa_w) - (mu+lambda_z)*z;
sys(2) = beta_y*x/(x+kappa_x)*kappa_z/(z+kappa_z) - (mu+lambda_y)*y; 




% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,state,input)


sys = state(2); % y


% end mdlOutputs