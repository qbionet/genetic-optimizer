%% Dynamics of the circuit with negative feedback
function dZdt = ODE_negativeFeedback(t,Z,par_proc,par_growth)

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y  = par_proc.lambda_y;
    F0        = par_proc.F0;
    theta_opt = par_proc.x_opt;
    sigma     = par_proc.sigma;  

    beta_w    = par_proc.beta_w;
    beta_z    = par_proc.beta_z;
    lambda_w  = par_proc.lambda_w;
    lambda_z  = par_proc.lambda_z;
    kappa_x   = par_proc.kappa_x;
    kappa_w   = par_proc.kappa_w;
    kappa_z   = par_proc.kappa_z;

    dZdt = zeros(3,1);

    w = Z(1);
    z = Z(2);
    y = Z(3);

    y_prod_rate = F0*exp(-(z-theta_opt)^2/(2*sigma^2));

    dZdt(1) = beta_w*x/(x+kappa_x)*kappa_z/(z+kappa_z) - (mu+lambda_w)*w;
    dZdt(2) = beta_z*w/(w+kappa_w) - (mu+lambda_z)*z;
    dZdt(3) = y_prod_rate - (mu+lambda_y)*y; 
end