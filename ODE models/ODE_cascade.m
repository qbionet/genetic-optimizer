%% Dynamics of a regulatory cascade
function dZdt = ODE_cascade(t,Z,par_proc,par_growth)

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y = par_proc.lambda_y;
    F0       = par_proc.F0;
    x_opt    = par_proc.x_opt;
    sigma    = par_proc.sigma;  

    beta_z   = par_proc.beta_z;
    lambda_z = par_proc.lambda_z;
    kappa_z  = par_proc.kappa_z;
    N        = par_proc.N;


    dZdt = zeros(N+1,1);

    for i = 1:N
        if i == 1
            dZdt(i) = beta_z*x/(x+kappa_z) - (mu + lambda_z)*Z(i);
        else
            dZdt(i) = beta_z*Z(i-1)/(Z(i-1)+kappa_z) - (mu + lambda_z)*Z(i);
        end
    end
    
    if N == 0
        y_prod_rate = F0*exp(-(x-x_opt)^2/(2*sigma^2));
    else
        y_prod_rate = F0*exp(-(Z_end-x_opt)^2/(2*sigma^2));
    end

    dZdt(N+1) = y_prod_rate - (mu+lambda_y)*y;      
end