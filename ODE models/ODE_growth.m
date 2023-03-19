%% Dynamics of growth regulation
function dZdt = ODE_growth(t,Z,par_proc,par_growth)

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y  = par_proc.lambda_y;

    beta_p    = par_proc.beta_p;
    beta_s    = par_proc.beta_s;
    beta_y    = par_proc.beta_y;
    lambda_p  = par_proc.lambda_p;
    lambda_s  = par_proc.lambda_s;
    kappa_x   = par_proc.kappa_x;
    eta_p     = par_proc.eta_p;
    eta_s     = par_proc.eta_s;
    omega     = par_proc.omega;

    dZdt = zeros(3,1);

    p = Z(1);
    s = Z(2);
    y = Z(3);

    dZdt(1) = beta_p - (mu+lambda_p+omega*s)*p;
    dZdt(2) = beta_s*x/(x+kappa_x) - (mu+lambda_s)*s; 
    dZdt(3) = beta_y*1/(1+p/eta_p)*1/(1+s/eta_s) - (mu+lambda_y)*y; 
end