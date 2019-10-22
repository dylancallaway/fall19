function [cost] = myPrinter1(theta_dot, delta_v_d)
t0 = tic;

% init variables defined in problem
T = .5;
dt = .001;
N_t = T/dt; % number of time steps
N_e = N_t; % number of droplets
g = -9.81;
v_2 = .25;
rho_1 = 2000;
rho_2 = 7000;
q_1 = 0;
q_2 = .001;
rho_star = (1-v_2)*rho_1 + v_2*rho_2;
q_star = (1-v_2)*q_1 + v_2*q_2;
R = .001;
V_e = (4/3)*pi*R^3;
m_e = V_e*rho_star;
q_e = V_e*q_star;
F_grav = [0, m_e*g, 0];
epsilon = 8.854e-12;
A_e = pi*R^2;
v_f = [.5, 0, .5];
rho_a = 1.225;
mu_f = 1.8e-5;
q_p = -8e-5;

time_steps = linspace(0, T, N_t);

theta_0 = [pi/2, 0, 0];
r_0 = [0, .5, 0]; % fixed arm end position (not the dispenser)
L = [.3, .2, .08];
L_bed = .8;

N_c = 10;
r_p_1 = linspace(-L_bed/2, L_bed/2, N_c)';
r_p(:, 3) = repmat(r_p_1, N_c, 1);
r_p_2 = [];
for i = 1:N_c
    r_p_2 = [r_p_2; repmat(r_p_1(i), N_c, 1)];
end
r_p(:, 1) = r_p_2;

last_stopped_drops = [0, 0, 0];

for n = 1:N_t
    % elapsed time
    t = time_steps(n);
    
    % dispenser state
    % integrate constant dispenser arm angular velocity to get angular
    % position
    theta = theta_0 + theta_dot*t;
    
    % _d is for dispenser
    % dispenser xyz position from dispenser angular position
    x_d = L(1)*cos(theta(1)) + L(2)*cos(theta(2)) + L(3)*sin(theta(3)); % pretty sure it should be L(3)*cos(theta(3)) and theta(2) is defined incorrectly in the diagram
    y_d = L(1)*sin(theta(1)) + L(2)*sin(theta(2));
    z_d = L(3)*cos(theta(3)); % should be L(3)*sin(theta(3))
    r_d = r_0 + [x_d, y_d, z_d];
    r_d_time(n, :) = r_d;
    
    % dispenser xyz velocity
    x_d_dot = -L(1)*theta_dot(1)*sin(theta(1)) - L(2)*theta_dot(2)*sin(theta(2)) + L(3)*theta_dot(3)*cos(theta(3));
    y_d_dot = L(1)*theta_dot(1)*cos(theta(1)) + L(2)*theta_dot(2)*cos(theta(2));
    z_d_dot = -L(3)*theta_dot(3)*sin(theta(3));
    v_d = [x_d_dot, y_d_dot, z_d_dot];
    
    % _e is for droplets
    % initialize new droplets
    r_e(n, :) = r_d;
    v_e(n, :) = v_d + delta_v_d;
    t_release(n, 1) = t;
    
    % droplet acceleration
    % electrical force
    r_e_temp = reshape(r_e, n, 1, 3);
    r_p_temp = reshape(r_p, 1, N_c^2, 3);
    dist = r_e_temp - r_p_temp;
    dist_mag = sqrt( dist(:,:,1).^2 + dist(:,:,2).^2 + dist(:,:,3).^2 );
    F_elec_temp = (q_p.*q_e)./(4.*pi.*epsilon.*dist_mag.^2).*dist; % droplet x point charge x dim
    F_elec_temp = sum(F_elec_temp, 2);
    F_elec = reshape(F_elec_temp, n, 3);
    
    % drag force
    Re = 2*R*rho_a*vecnorm(v_f-v_e, 2, 2)/mu_f;
    C_D(Re>0 & Re<=1) = 24./Re(Re>0 & Re<=1);
    C_D(Re>1 & Re<=400) = 24./(Re(Re>1 & Re<=400).^.0646);
    C_D(Re>400 & Re<=3e5) = .5;
    C_D(Re>3e5 & Re<=2e6) = .000366.*(Re(Re>3e5 & Re<=2e6).^.4275);
    C_D(Re>2e6) = .18;
    
    F_drag = .5.*rho_a.*C_D'.*vecnorm(v_f - v_e, 2, 2).*(v_f - v_e).*A_e;
    
    % total force
    F_e = F_grav + F_elec + F_drag;
    a_e = F_e/m_e;
    
    stopped_drops = repmat(r_e(:,2)<=0, 1, 3);
    r_e(~stopped_drops) = r_e(~stopped_drops) + dt*v_e(~stopped_drops);
    v_e(~stopped_drops) = v_e(~stopped_drops) + dt*a_e(~stopped_drops);
    v_e(stopped_drops) = 0;
    a_e(stopped_drops) = 0;
    
    land_check = ~(stopped_drops == last_stopped_drops);
    t_land(n,1) = 0;
    t_land(land_check(:,1)) = t;
    
    last_stopped_drops = [stopped_drops; 0, 0, 0];
    
    %     disp(n)
    t1 = toc(t0)
end
r_des = rand(N_t, 3);
numer = sum(vecnorm(r_des - r_e, 2, 2), 1);
denom = sum(vecnorm(diff(r_des, 1, 1), 2, 2), 1);
cost = numer/denom;

if 1
    % part 1
    % final pattern plot (number 2)
    figure;
    tof = t_land - t_release;
    landed = tof>0;
    r_e_x = r_e(:,1);
    r_e_z = r_e(:,3);
    landed_scat = scatter(r_e_x(landed), r_e_z(landed), 100, tof(landed), '.');
    hold on
    grid on
    c = colorbar;
    c.Label.String = 'Time of Flight of Landed Droplets (s)';
    not_landed_scat = scatter(r_e_x(~landed), r_e_z(~landed), 1, 'r','.');
    grid_scat = scatter(r_p(:,1), r_p(:,3), 5, 'm','o');
    xlabel('X Position (m)')
    ylabel('Z Position (m)')
    legend([not_landed_scat, grid_scat], 'Droplets in Flight', 'Charges on Substrate')
    hold off
    
    % part 2
    figure;
    r_e_y = r_e(:,2);
    r_e_y(landed) = 0;
    dispenser_scat = scatter3(r_d_time(:,1), r_d_time(:,3), r_d_time(:,2), 10, 'k.');
    hold on
    xlabel('X Position (m)')
    ylabel('Z Position (m)')
    zlabel('Y Position (m)')
    dispenser_start_scat = scatter3(r_d_time(1,1), r_d_time(1,3), r_d_time(1,2), 50, 'gx');
    dispenser_end_scat = scatter3(r_d_time(end,1), r_d_time(end,3), r_d_time(end,2), 50, 'cx');

    zlim([0, max(r_d_time(:,2))])
    landed_scat = scatter(r_e_x(landed), r_e_z(landed), 100, tof(landed), '.');
    not_landed_scat = scatter3(r_e_x(~landed), r_e_z(~landed), r_e_y(~landed), 1, 'r','.');
    grid_scat = scatter(r_p(:,1), r_p(:,3), 5, 'm','o');
    legend([dispenser_scat, not_landed_scat, grid_scat, dispenser_start_scat, dispenser_end_scat], 'Dispenser Path', 'Droplets in Flight', 'Charges on Substrate', 'Dispenser Start Point', 'Dispenser End Point')
    hold off
    
    % part 3
    % part 4
    % just writing and setting electircal force to 0
end
end