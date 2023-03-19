%% Dynamics of a single reporter
function dZdt = ODE_reporter(t,Z,par_proc,par_growth)

    % Process to be regulated (this needs to be changed when considering a different example)
    y   = Z(1);

    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;

    % Parameters of the process (these need to be defined in process_parameters.m)
    lambda_y = par_proc.lambda_y;
    F0       = par_proc.F0;
    x_opt    = par_proc.x_opt;
    sigma    = par_proc.sigma; 


    dZdt = zeros(1,1);

    dZdt(1) = F0*exp(-(x-x_opt)^2/(2*sigma^2)) - (mu+lambda_y)*y;
       
end