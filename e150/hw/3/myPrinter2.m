function [cost] = myPrinter2(theta_dot, delta_v_d, plots)
t0 = tic;

% init variables defined in problem
T = 3;
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
t = time_steps';

theta_0 = [pi/2, 0, 0];
r_0_y = [0, .5, 0]; % fixed arm end position (not the dispenser)
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


% dispenser state
% integrate constant dispenser arm angular velocity to get angular
% position
theta = theta_0 + theta_dot.*t;

% _d is for dispenser
% dispenser xyz position from dispenser angular position
x_d = L(1).*cos(theta(:,1)) + L(2).*cos(theta(:,2)) + L(3).*sin(theta(:,3)); % pretty sure it should be L(3)*cos(theta(3)) and theta(2) is defined incorrectly in the diagram
y_d = L(1).*sin(theta(:,1)) + L(2).*sin(theta(:,2));
z_d = L(3).*cos(theta(:,3)); % should be L(3)*sin(theta(3))
r_d = r_0_y + [x_d, y_d, z_d]; % row is time step x col is dim

% dispenser xyz velocity
x_d_dot = -L(1).*theta_dot(1).*sin(theta(:,1)) - L(2).*theta_dot(2).*sin(theta(:,2)) + L(3).*theta_dot(3).*cos(theta(:,3));
y_d_dot = L(1).*theta_dot(1).*cos(theta(:,1)) + L(2).*theta_dot(2).*cos(theta(:,2));
z_d_dot = -L(3).*theta_dot(3).*sin(theta(:,3));
v_d = [x_d_dot, y_d_dot, z_d_dot]; % row is time step x col is dim

% _e is for droplets
% initialize new droplets
r_e_release = r_d;
v_e_release = v_d + delta_v_d;
t_release = t;

% droplet acceleration
% total force
F_e = repmat(F_grav, N_t, 1);
a_e = F_e/m_e;

r_0_y = r_e_release(:,2);
v_0_y = v_e_release(:,2);
a_e_y = a_e(:,2);

tof = (-v_0_y - sqrt( v_0_y.^2 - 2.*a_e_y.*r_0_y )) ./ (a_e_y);

% t_land = t_release + tof;
% landed = t_land<=T;
% tof(~landed) = max(t_land) - t_land(~landed);

r_e_final = r_e_release + v_e_release.*tof + .5.*a_e.*tof.^2;

% r_des = rand(N_t, 3);
r_des_raw = load("robotprint_data.mat");
r_des = r_des_raw.ri(:,1:N_t)';
numer = sum(vecnorm(r_des - r_e_final, 2, 2), 1);
denom = sum(vecnorm(diff(r_des, 1, 1), 2, 2), 1);
cost = numer/denom;

if plots
    % part 3
    % final pattern plot (number 2)
    figure;
%     landed = t_land<=T;
    r_e_x = r_e_final(:,1);
    r_e_z = r_e_final(:,3);
    r_des_x = r_des(:,1);
    r_des_z = r_des(:,3);
    des_scat = scatter(r_des_x, r_des_z, 25, 'go', 'filled');
    hold on
    landed_scat = scatter(r_e_x, r_e_z, 5, 'bo', 'filled');
    grid on
%     c = colorbar;
%     c.Label.String = 'Time of Flight of Droplets (s)';
    %     grid_scat = scatter(r_p(:,1), r_p(:,3), 5, 'm','o');

    xlabel('X Position (m)')
    ylabel('Z Position (m)')
    legend('Desired Droplet Pattern', 'Generated Droplet Pattern')
    hold off
    
    % part none :/
%     figure;
%     r_e_y = r_e_final(:,2);
%     r_e_y = 0;
%     dispenser_scat = scatter3(r_d(:,1), r_d(:,3), r_d(:,2), 10, 'k.');
%     hold on
%     xlabel('X Position (m)')
%     ylabel('Z Position (m)')
%     zlabel('Y Position (m)')
%     dispenser_start_scat = scatter3(r_d(1,1), r_d(1,3), r_d(1,2), 50, 'gx');
%     dispenser_end_scat = scatter3(r_d(end,1), r_d(end,3), r_d(end,2), 50, 'cx');
%     zlim( [0, max( [max(r_e_y), max(r_d(:,2))] )] )
%     landed_scat = scatter(r_e_x, r_e_z, 100, tof, '.');
%     r_des_y = r_des(:,2);
%     des_scat = scatter3(r_des_x, r_des_z, r_des_y, 5, 'mo');
%     legend([dispenser_scat, des_scat, dispenser_start_scat, dispenser_end_scat], 'Dispenser Path', 'Desired Droplet Pattern', 'Dispenser Start Point', 'Dispenser End Point')
%     hold off
end
end