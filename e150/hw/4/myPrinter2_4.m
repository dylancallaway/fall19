function [test] = myPrinter2_4(plots)
% variable glossary for part 1
% mu_fo = 0.001; % viscosity of fluid
sigma_f0 = 0.617;
sigma_p0 = 0.13;
% rho_f0 = 2000;
% rho_p0 = 4000;
C_f = 1600;
C_p = 3800;
% c_1 = 0.01;
% c_2 = 2;
% R = 0.001;
% Q_0 = 1e-6;
% b = 0.25;
% omega = 12;
h = 10;
% theta_0 = 300;
theta_a = 300;
% theta_dot = 50;
% dt = 0.01;
% T = 2.5;
% a = 0.8;
% k1 = 0.5;
% k2 = 1;
% k3 = 2;
% v_2_all = [0, 0.05, 0.1, 0.15, 0.2, 0.25]; % vol frac of particles
phi = 0.5;

% part 2
R_0 = 0.001;
q = 2;
v_m0 = 0.01; % initial mean velocity
Q_0 = pi*R_0^2*v_m0;
mu_f0 = 0.001;
rho_f0 = 2000;
rho_p0 = 4000;
theta_bar = 500;
v_2_all = [0, 0.05, 0.1, 0.15, 0.2, 0.25]; % vol frac of particles
theta_0 = 300;
tau_crit = 0.001;
eta = (10^-5)/3600;
c_1 = 0.01;
c_2 = 2;
dt = 3.6;
T = 3600*100;
k1 = 0.5;
k2 = 1;
k3 = 2;

theta = theta_bar;

% desired flow rate
Q = Q_0;

% R_ss
R_ss_norm_all = [];

ind_1 = 1;
p = [];
v_2 = .15;
theta = 100*(3:8);



norm_R_ss_all = [];
v_2_all = [0, 0.05, 0.1, 0.15, 0.2, 0.25]; % vol frac of particles

for v_2 = v_2_all
    % other inits
    t = 0;
    t_all = [];
    dp_dx_all = [];
    Re_all = [];
    R = R_0;
    
    norm_R_all = [];
    
    v_p = v_2;
    v_f = 1 - v_p;
    
    % viscosity
    mu_f = mu_f0.*exp(-k1.*(theta-theta_0)./theta_0);
    
    % conductivity of fluid
    sigma_f = sigma_f0.*exp(-k2.*(theta-theta_0)./theta_0);
    
    % conductivity of particles
    sigma_p = sigma_p0.*exp(-k3.*(theta-theta_0)./theta_0);
    
    % effective viscosity
    mu_star = mu_f.*(1+2.5.*v_p./(1-v_p));
    
    % effective density
    rho_star = v_f.*rho_f0 + v_p.*rho_p0;
    
    % R_ss_all
    R_ss = ( mu_star.*Q.*(q+2)./pi./tau_crit ).^(1/3);
    R_ss_norm_theta = R_ss./R_0;
    
    % effective heat capacity
    C_star = v_f.*C_f + v_p.*C_p;
    
    %     norm_R_ss_all(ind_1, :) = R_ss_norm_theta
    while t <= T
        % gamma_star
        gamma_star = 2.*c_1.*Q.*rho_star./pi./R./mu_star;
        
        % velocity profile characteristics
        %         q_plus = 0.5.*( (gamma_star+c_2) + sqrt( (gamma_star + c_2).^2 + 8.*gamma_star ) );
        %         q_minus = 0.5.*( (gamma_star+c_2) - sqrt( (gamma_star + c_2).^2 + 8.*gamma_star ) );
        %         q = q_plus;
        v_max = Q.*(q+2)./pi./R.^2./q;
        Re = rho_star.*v_max.*2.*R./mu_star;
        Re_all = [Re_all, Re];
        
        % wall shear stress
        tau_w = mu_star.*Q.*(q+2)./pi./R.^3;
        
        % pressure gradient
        C = 2.*mu_star.*(q+2)./pi./R.^4;
        dp_dx = C.*Q;
        dp_dx_all = [dp_dx_all, dp_dx];
        
        % effective conductivity
        sigma_star_plus = sigma_p + (1-v_p)./( (1./(sigma_f-sigma_p)) + (v_p./3./sigma_p) );
        sigma_star_minus = sigma_f + (v_p)./((1./(sigma_p-sigma_f)) + (1-v_p)./3./sigma_f);
        sigma_star = phi.*sigma_star_minus + (1-phi).*sigma_star_plus;
        
        % joule heating
        next_theta = theta;
        S = h.*(theta-theta_a)./R;
        % current
        %         J = sqrt( (sigma_star/a).*( (rho_star.*C_star.*(next_theta-theta)./dt) + S ) );
        %         J_all = [J_all, J];
        
        % nozzle radius
        dR_dt = eta.*max( (mu_star.*Q.*(q+2)./pi./tau_crit./R.^3) - 1 , 0);
        next_R = R + dt.*dR_dt;
        norm_R = R./R_0;
        
        % increment time
        t = t + dt;
        t_all = [t_all, t];
        R = next_R;
    end
    
    %     norm_R_all = [norm_R_all, norm_R];
    
    if plots
        colors = lines(6);
        % part 2, 4
        figure(1)
        p(ind_1) = plot(theta, R_ss_norm_theta, '-', 'Color', colors(ind_1,:));
        hold on
%         p2(ind_1) = plot(theta, norm_R, '--', 'Color', colors(ind_1,:));
        %         yline(R_ss_norm, '--', 'Color', colors(ind_1,:));
        xlabel('Temperature (K)')
        ylabel('Normalized Radius (R_{ss}/R_0)')
        if ind_1 == 6
            legend([p(1), p(2), p(3), p(4), p(5), p(6)], 'v_2 = 0', 'v_2 = 0.05', 'v_2 = 0.1', 'v_2 = 0.15', 'v_2 = 0.2', 'v_2 = 0.25')
        end
        
        % part 2, number 5
        % percent difference for this theta
        if ind_1 == 1
        fprintf("Temperature & 300 & 400 & 500 & 600 & 700 & 800\n")
        end
        fprintf("v_2 & %.2f", v_2)
        for p_ind = 1:6
           pct_diff = 100*(R_ss_norm_theta(p_ind)-norm_R(p_ind))/R_ss_norm_theta(p_ind);
           fprintf(" & %.2f", pct_diff) 
        end
        fprintf("\n")
        
    end
    
    ind_1 = ind_1 + 1;
end
% figure(2)
% contourf(theta, v_2_all, norm_R_all, 100);
% colorbar()
test = 0;
end

