function [sys,x0,str,ts] = reporter_Sfun(t,state,input,flag)



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
sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);


x0  = rand()*100;
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
F0       = par_proc.F0;          % height of objective function
sigma    = par_proc.sigma;       % width of objective function
lambda_y = par_proc.lambda_y;    % degradation rate constant of the reporter (1/h)

% Parameters that are shared: growth rate
mu       = par_growth.mu;       % growth rate (1/h)

y     = state;
x     = input(1);
x_opt = input(2);


% Process (this needs to be changed when considering a different example)
sys(1)  = F0*exp(-(x-x_opt)^2/(2*sigma^2)) - (mu+lambda_y)*y;



% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,state,input)

y = state;

sys = y;


% end mdlOutputs