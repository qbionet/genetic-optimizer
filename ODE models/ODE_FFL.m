%% Dynamics of the feedforward motif
function dZdt = ODE_FFL(t,Z,par_proc,par_growth)

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y  = par_proc.lambda_y;

    beta_z    = par_proc.beta_z;
    beta_y    = par_proc.beta_y;
    lambda_z  = par_proc.lambda_z;
    kappa_x   = par_proc.kappa_x;
    kappa_z   = par_proc.kappa_z;
    kappa_w   = par_proc.kappa_w;
    w         = par_proc.w;         % treated as input in the Simulink realization

    dZdt = zeros(2,1);

    z = Z(1);
    y = Z(2);

    dZdt(1) = beta_z*(x/kappa_x)/(1 + w/kappa_w) - (mu+lambda_z)*z;
    dZdt(2) = beta_y*x/(x+kappa_x)*kappa_z/(z+kappa_z) - (mu+lambda_y)*y; 
end