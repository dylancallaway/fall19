dbstop if error
format compact

%% Problem 1
close all
clear
clc

% a
C_D_p = .0356;
AR = 8.75;
S = .69;
m = 5.3;
g = 9.81;
f_L = m.*g;
rho  = 1.225;
v = 3:.1:30;

f_D_p = .5.*rho.*v.^2.*S.*C_D_p;
C_L = 2.*f_L./rho./v.^2./S;
C_D_i = C_L.^2/pi/AR;
f_D_i = .5.*rho.*v.^2.*S.*C_D_i;
f_D_t = f_D_p + f_D_i;
P_t = f_D_t .* v;

% b
v_feas_log = f_D_t<=10 & P_t<=200;
v_feas = v(v_feas_log);
v_feas_bound = [min(v_feas), max(v_feas)] % m/s

% c
v_r_ind = f_D_t == min(f_D_t);
v_r = v(v_r_ind);

v_e_ind = P_t == min(P_t);
v_e = v(v_e_ind);

% a
figure
subplot(1, 2, 1)
plot(v, f_D_p)
hold on
plot(v, f_D_i)
plot(v, f_D_t)
a1 = area(v_feas, f_D_t(v_feas_log));
a1.FaceAlpha = .1;
a1.EdgeAlpha = 0;
a1.FaceColor = 'g';
plot(v_r, f_D_t(v_r_ind), 'blackx')
hold off
legend('f_D_p (Parasitic)', 'f_D_i (Induced)', 'f_D_t (Total)', 'Feasible Air Speeds', 'Maximum Range Speed')
xlabel('Air Speed (m/s)')
ylabel('Drag Force (N)')
title('Drag Force vs. Air Speed')

subplot(1, 2, 2)
plot(v, P_t)
hold on
a1 = area(v_feas, P_t(v_feas_log));
a1.FaceAlpha = .1;
a1.EdgeAlpha = 0;
a1.FaceColor = 'g';
plot(v_e, P_t(v_e_ind), 'blackx')
hold off
legend('Total Drag Power', 'Feasible Air Speeds', 'Maximum Endurance Speed')
xlabel('Air Speed (m/s)')
ylabel('Drag Power (W)')
title('Drag Power vs. Air Speed')

%% Problem 2
%% Problem 3
%% Problem 4
close all
clear
clc

alpha_s = [5; 0; 10]
T_BS = [sqrt(2)/2, sqrt(2)/2, 0; -sqrt(2)/2, sqrt(2)/2, 0; 0, 0, 1]
alpha_B = T_BS * alpha_s % m/s^2

%% Problem 5
%% Problem 6


