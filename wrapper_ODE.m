clear all
close all
clc


% load parameters of the optimizer
load("Parameters/optimizer_parameters.mat")
load("Parameters/process_parameters.mat")
load("Parameters/growth_parameters.mat")




T_end = 100;                % duration of simulation
x0 = rand(1,22 + par_proc.dim)*100;        % random initial conditions (dimension is 22+process)
[tz,z] = ode45(@(t,z) ODE_closedLoop(t,z,par_opt,par_proc,par_growth), [0 T_end], x0);

figure()
hold on
plot([0 T_end],par_proc.x_opt*ones(1,2),'r','LineWidth',5)
plot(tz,z(:,21),'b','LineWidth',5)
ylim([0 150]), grid on
legend('optimum','closed loop')
xlabel('time (h)'), ylabel('x (nM)')


function dzdt = ODE_closedLoop(t,z,par_opt,par_proc,par_growth)

    dzdt = zeros(22 + par_proc.dim,1);

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

    % Process to be regulated 
    % Extra variables have to be added for custom examples
    y   = z(23);


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

    % Parameters of the process (this needs to be changed when considering a different example)
    lambda_y = par_proc.lambda_y;
    F0       = par_proc.F0;
    x_opt    = par_proc.x_opt;
    sigma    = par_proc.sigma; 

    % Parameters shared: growth
    mu       = par_growth.mu;



    % Delay    
    dzdt(1)  = beta_xd*x/(x+K_x) - (mu + lambda_d)*xd;
    dzdt(2)  = beta_yd*y/(y+K_y) - (mu + lambda_d)*yd;

    % Repressilator
    dzdt(3)  = beta_o/(1+(r/K_o)^4) - (mu+lambda_o)*o1;
    dzdt(4)  = beta_o/(1+(o1/K_o)^2) - (mu+lambda_o)*o2;
    dzdt(5)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*r;
    dzdt(6)  = beta_o/(1+(o2/K_o)^2) - (mu+lambda_o)*a;

    % gRNAs
    dzdt(7)  = alpha_x*K_c^4/(r^4+K_c^4)*xd/(xd+K_xb) + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xm;
    dzdt(8)  = alpha_x*K_c^4/(r^4+K_c^4)*x/(x+K_xb)   + alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xp;
    dzdt(9)  = alpha_y*K_c^4/(r^4+K_c^4)*yd/(yd+K_yb) + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ym;
    dzdt(10) = alpha_y*K_c^4/(r^4+K_c^4)*y/(y+K_yb)   + alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*yp;

    % STARs
    dzdt(11)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xp^2 + k_xb^2) - (mu + delta_c)*xmb;
    dzdt(12)  = alpha_xb*a^4/(a^4+K_c^4)*k_xb^2/(xm^2 + k_xb^2) - (mu + delta_c)*xpb;
    dzdt(13)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(yp^2 + k_yb^2) - (mu + delta_c)*ymb;
    dzdt(14)  = alpha_yb*a^4/(a^4+K_c^4)*k_yb^2/(ym^2 + k_yb^2) - (mu + delta_c)*ypb;

    % AND gates
    dzdt(15) = alpha_q*xpb/(xpb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qpp;
    dzdt(16) = alpha_q*xmb/(xmb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qmm;
    dzdt(17) = alpha_q*xpb/(xpb+k_xt)*ymb/(ymb+k_yt) - (mu + delta_q)*qpm;
    dzdt(18) = alpha_q*xmb/(xmb+k_xt)*ypb/(ypb+k_yt) - (mu + delta_q)*qmp;

    % OR gates
    dzdt(19) = beta_u*(qpp/(k_q + qpp)+qmm/(k_q + qmm)) - (mu + lambda_u)*u1;
    dzdt(20) = beta_u*(qpm/(k_q + qpm)+qmp/(k_q + qmp)) - (mu + lambda_u)*u2;

    % Regulator
    dzdt(21) = beta_x*x/(x+K_x) + beta_xb*u1^2/(u1^2+K_u1^2) - (mu + lambda_x*(1+nu*v))*x;
    
    % Protease
    dzdt(22) = beta_v*u2/(u2+K_u2) - (mu + lambda_v)*v;

    % Process (this needs to be changed for custom examples) 
    % ODE models for examples featured in the paper can be found in "ODE models")
    % The corresponding parameters need to be added to
    % "Parameters/process_parameters.mat" and then run so that they are
    % saved and can be loaded to the workspace
    dzdt(23) = F0*exp(-(x-x_opt)^2/(2*sigma^2)) - (mu+lambda_y)*y;
       
end
