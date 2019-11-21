function [test] = myPrinter1(plots)
% variable glossary for part 1
mu_f0 = 0.001; % viscosity of fluid
sigma_f0 = 0.617;
sigma_p0 = 0.13;
rho_f0 = 2000;
rho_p0 = 4000;
C_f = 1600;
C_p = 3800;
c_1 = 0.01;
c_2 = 2;
R = 0.001;
Q_0 = 1e-6;
b = 0.25;
omega = 12;
h = 10;
theta_0 = 300;
theta_a = 300;
theta_dot = 50;
dt = 0.01;
T = 2.5;
a = 0.8;
k1 = 0.5;
k2 = 1;
k3 = 2;
v_2_all = [0, 0.05, 0.1, 0.15, 0.2, 0.25]; % vol frac of particles
phi = 0.5;



for v_2 = v_2_all
    % other inits
    theta = theta_0;
    I_all = [];
    t = 0;
    t_all = [];
    dp_dx_all = [];
    Re_all = [];
    
    v_p = v_2;
    v_f = 1 - v_p;
    while t <= T
        % viscosity at temperature
        mu_f = mu_f0.*exp(-k1.*(theta-theta_0)./theta_0);
        
        % conductivity of fluid at temp
        sigma_f = sigma_f0.*exp(-k2.*(theta-theta_0)./theta_0);
        
        % conductivity of particles at temp
        sigma_p = sigma_p0.*exp(-k3.*(theta-theta_0)./theta_0);
        
        % effective viscosity
        mu_star = mu_f.*(1+2.5.*v_p./(1-v_p));
        
        % effective density
        rho_star = v_f.*rho_f0 + v_p.*rho_p0;
        
        % desired flow rate
        Q = Q_0.*(1+b.*sin(omega.*t));
        
        % gamma_star
        gamma_star = 2.*c_1.*Q.*rho_star/pi/R/mu_star;
        
        % velocity profile characteristics
        q_plus = 0.5.*( (gamma_star+c_2) + sqrt( (gamma_star + c_2).^2 + 8.*gamma_star ) );
        q_minus = 0.5.*( (gamma_star+c_2) - sqrt( (gamma_star + c_2).^2 + 8.*gamma_star ) );
        q = q_plus;
        v_max = Q*(q+2)/pi/R^2/q;
        Re = rho_star*v_max*2*R/mu_star;
        Re_all = [Re_all, Re];
        
        % wall shear stress
        tau_w = mu_star.*Q.*(q+2)/pi/R.^3;
        
        % pressure gradient
        C = 2.*mu_star.*(q+2)/pi/R.^4;
        dp_dx = -C.*Q;
        dp_dx_all = [dp_dx_all, dp_dx];
        
        % effective heat capacity
        C_star = v_f.*C_f + v_p.*C_p;
        
        % effective conductivity
        sigma_star_plus = sigma_p + (1-v_p)./( (1./(sigma_f-sigma_p)) + (v_p/3/sigma_p) );
        sigma_star_minus = sigma_f + (v_p)./((1/(sigma_p-sigma_f)) + (1-v_p)/3/sigma_f);
        sigma_star = phi.*sigma_star_minus + (1-phi).*sigma_star_plus;
        
        % joule heating
        next_theta = theta + theta_dot*dt;
        S = h.*(theta-theta_a)./R;
        % current
        J = sqrt( (sigma_star/a).*( (rho_star.*C_star.*(next_theta-theta)./dt) + S ) );
        I = J*pi*R^2;
        I_all = [I_all, I];
        
        % increment time
        t = t + dt;
        t_all = [t_all, t];
        theta = next_theta;
    end
    if plots
        % part 1, 6
        figure(1)
        plot(t_all, I_all)
        xlabel('Time (s)')
        ylabel('Current Through Mixture (A)')
        hold on
        legend('v_2 = 0', 'v_2 = 0.05', 'v_2 = 0.1', 'v_2 = 0.15', 'v_2 = 0.2', 'v_2 = 0.25')
        
        % part 1, 7
        figure(2)
        plot(t_all, dp_dx_all)
        xlabel('Time (s)')
        ylabel('Pressure Gradient Along Nozzle (Pa/m)')
        hold on
        legend('v_2 = 0', 'v_2 = 0.05', 'v_2 = 0.1', 'v_2 = 0.15', 'v_2 = 0.2', 'v_2 = 0.25')
        
        % part 1, 8
        figure(3)
        plot(t_all, Re_all)
        xlabel('Time (s)')
        ylabel('Reynolds Number of Flow')
        hold on
        legend('v_2 = 0', 'v_2 = 0.05', 'v_2 = 0.1', 'v_2 = 0.15', 'v_2 = 0.2', 'v_2 = 0.25')
    end
end
test = 0;
end

