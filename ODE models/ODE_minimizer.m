%% Dynamics of minimizer
function dZdt = ODE_minimizer(t,Z,par_proc,par_growth)

    % y needs to be provided by the dynamics of the process

    
    % Parameters shared: growth (this needs to be defined in growth_parameters.m)
    mu       = par_growth.mu;


    % Optimizer
    xd  = Z(1);
    yd  = Z(2);
    o1  = Z(3);
    o2  = Z(4);
    r   = Z(5);
    a   = Z(6);
    xm  = Z(7);  % gRNA
    xp  = Z(8);  % gRNA
    ym  = Z(9);  % gRNA
    yp  = Z(10); % gRNA
    xmb = Z(11); % STAR
    xpb = Z(12); % STAR
    ymb = Z(13); % STAR
    ypb = Z(14); % STAR
    qpp = Z(15); 
    qmm = Z(16);
    qpm = Z(17);
    qmp = Z(18);
    u1  = Z(19);
    u2  = Z(20);
    x   = Z(21);
    v   = Z(22);

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

    k_x      = par_opt.k_x;
    k_y      = par_opt.k_y;
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

    dZdt = zeros(3,1);

    % Delay    
    dZdt(1)  = beta_xd*x/(x+K_x) - (mu + lambda_d)*xd;
    dZdt(2)  = beta_yd*y/(y+K_y) - (mu + lambda_d)*yd;

    % Repressilator
    dZdt(3)  = beta_o/(1+(r/K_o)^4) - (mu+lambda_o)*o1;
    dZdt(4)  = beta_o/(1+(o1/K_o)^2) - (mu+lambda_o)*o2;
    dZdt(5)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*r;
    dZdt(6)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*a;

    % gRNAs
    dZdt(7)  = alpha_x*K_c^4/(r^4+K_c^4)*xd/(xd+K_xb) + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xm;
    dZdt(8)  = alpha_x*K_c^4/(r^4+K_c^4)*x/(x+K_xb)   + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xp;
    dZdt(9)  = alpha_y*K_c^4/(r^4+K_c^4)*yd/(yd+K_yb) + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ym;
    dZdt(10) = alpha_y*K_c^4/(r^4+K_c^4)*y/(y+K_yb)   + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*yp;

    % STARs
    dZdt(11)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xmb;
    dZdt(12)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xpb;
    dZdt(13)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ymb;
    dZdt(14)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*ypb;

    % AND gates
    dZdt(15) = alpha_q*xpb/(xpb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qpp;
    dZdt(16) = alpha_q*xmb/(xmb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qmm;
    dZdt(17) = alpha_q*xpb/(xpb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qpm;
    dZdt(18) = alpha_q*xmb/(xmb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qmp;

    % OR gates
    dZdt(19) = beta_u*(qpm/(k_q + qpm)+qmp/(k_q + qmp)) - (mu + lambda_u)*u1;
    dZdt(20) = beta_u*(qpp/(k_q + qpp)+qmm/(k_q + qmm)) - (mu + lambda_u)*u2;

    % Regulator
    dZdt(21) = beta_x*x/(x+K_x) + beta_xb*u1^2/(u1^2+K_u1^2) - (mu + lambda_x*(1+nu*v))*x;
    
    % Protease
    dZdt(22) = beta_v*u2/(u2+K_u2) - (mu + lambda_v)*v; 
end