%% Dynamics of a self-activated cascade
function dZdt = ODE_selfActivated(t,Z,par_proc,par_growth)

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y  = par_proc.lambda_y;
    F0        = par_proc.F0;
    theta_opt = par_proc.x_opt;
    sigma     = par_proc.sigma;  

    beta_z1   = par_proc.beta_z1;
    beta_z2   = par_proc.beta_z2;
    lambda_z  = par_proc.lambda_z;
    kappa_x   = par_proc.kappa_x;
    kappa_z   = par_proc.kappa_z;

    dZdt = zeros(2,1);

    z = Z(1);
    y = Z(2);

    y_prod_rate = F0*exp(-(z-theta_opt)^2/(2*sigma^2));

    dZdt(1) = beta_z1*x/(x+kappa_x) + beta_z2*z/(z+kappa_z) - (mu+lambda_z)*z;
    dZdt(2) = y_prod_rate - (mu+lambda_y)*y; 
end