%% HW3
% Dylan Callaway
% Fall 19 - ME102B

dbstop if error
format compact
close all
clear
clc

%% 1
r0 = 168/4;
r1 = 168*3/4;

d0 = r0*2 % mm
d1 = r1*2 % mm

m = 4; % mm

N0 = d0/m
N1 = d1/m

%% 2
clear
clc

N0 = 20
N1 = N0*2

P = 8;

d0 = N0/P; % in
d1 = N1/P; % in

center_to_center = (d0+d1)/2 % in

p = pi/P

%% 3
clear
clc

%% a
rpm_a = 1200; % rpm
N_a = 15;
N_b1 = 45;
N_b2 = 15;
N_c = 45;
rpm_b = rpm_a*N_a/N_b1 % rpm
rpm_c = rpm_b*N_b2/N_c % rpm

P = 5;
d_a = N_a/P % in
d_b1 = N_b1/P
d_b2 = N_b2/P
d_c = N_c/P

p = pi/P % ??

%% b
% Assume output shaft has torque? Otherwise forces are 0 for 100%
% efficiency? AKA assume operating at 1kW and 1200RPM, not rated values.
%% i
w_a = rpm_a/60*2*pi; % rad/s
w_b = rpm_b/60*2*pi;
w_c = rpm_c/60*2*pi;

P_a = 1e3; % W

T_a = P_a/w_a % Nm
T_b = P_a/w_b
T_c = P_a/w_c

%% ii
T_a_ii = T_a % Nm
P_b_ii = .95*P_a; % W
T_b_ii = P_b_ii/w_b % Nm
P_c_ii = .95*P_b_ii; % W
T_c_ii = P_c_ii/w_c % Nm

%% c
r_a = d_a/2/1000*25.4; % m
r_b1 = d_b1/2/1000*25.4;
r_b2 = d_b2/2/1000*25.4;
r_c = d_c/2/1000*25.4;

phi = 25; % deg
F_ab1_mag = T_a/r_a; % N
F_ab1 = F_ab1_mag.*[cosd(phi), -sind(phi), 0];
F_b2c_mag = T_c/r_c;
F_b2c = F_b2c_mag.*[-cosd(phi), -sind(phi), 0];
r_b1A = [0, 0, .025]; % m
r_BA = [0, 0, -.1];
r_cA = [0, 0, -.125];

syms A_x A_y A_z B_x B_y B_z

A = [A_x, A_y, A_z];
B = [B_x, B_y, B_z];
sum_F = A + B + F_ab1 + F_b2c == [0, 0, 0];
sum_M = cross(r_b1A, F_ab1) + cross(r_BA, B) + cross(r_cA, F_b2c) == [0, 0, 0];
[A_x, A_y, A_z, B_x, B_y, B_z] = solve([sum_F, sum_M], [A_x A_y A_z B_x B_y B_z]);

A = double([A_x, A_y, A_z]) % N
B = double([B_x, B_y, B_z])

%% 4
clear
clc

%% a
N_1 = 16;
N_2 = 32;
N_3 = 24;

T_1 = 100; % lbfin
T_2 = T_1*N_2/N_1;
T_3 = T_2*N_3/N_2;

P = 8;
phi = 20;

r_1 = N_1/P/2; % in
r_2 = N_2/P/2;
r_3 = N_3/P/2;

F_21_mag = T_1/r_1; % lbf
F_21 = F_21_mag.*[-cosd(phi), -sind(phi), 0]';

F_23_mag = T_3/r_3;
F_23 = F_23_mag.*[-sind(phi), -cosd(phi), 0]';

syms A_x A_y A_z B_x B_y B_z
A = [A_x; A_y; A_z];
B = [B_x; B_y; B_z];

r_BA = [0, 0, 2];
r = [0, 0, -1];

sum_F = F_21 + F_23 + A + B == 0;
sum_M = cross(r_BA, B) + cross(r, F_21) + cross(r, F_23) == 0;

[A_x, A_y, A_z, B_x, B_y, B_z] = solve([sum_F, sum_M], [A_x A_y A_z B_x B_y B_z]);

A = double([A_x; A_y; A_z]) % lbf
B = double([B_x; B_y; B_z])

%% b
F_21 = F_21_mag.*[cosd(phi), -sind(phi), 0]';

F_23 = F_23_mag.*[-sind(phi), cosd(phi), 0]';

syms A_x A_y A_z B_x B_y B_z
A = [A_x; A_y; A_z];
B = [B_x; B_y; B_z];

r_BA = [0, 0, 2];
r = [0, 0, -1];

sum_F = F_21 + F_23 + A + B == 0;
sum_M = cross(r_BA, B) + cross(r, F_21) + cross(r, F_23) == 0;

[A_x, A_y, A_z, B_x, B_y, B_z] = solve([sum_F, sum_M], [A_x A_y A_z B_x B_y B_z]);

A = double([A_x; A_y; A_z]) % lbf
B = double([B_x; B_y; B_z])

% The overall magnitude of the reaction forces are smaller because some of
% the loading that was supported by the bearings in part (a) is cancelled
% out by gear 3 reaction forces in part (b). This is good to note thatr the
% phyiscal arragnement you put your gears in may depend on what direction
% they typically spin.

%% 5
clear
clc

d_t = .025; % m
p = .005;
F = 5e3; % N
mu_c = .06;
mu_t = .09;
d_c = .045; % m
d_m = d_t - p/2



T_R = (F*d_m/2)*(p + pi*mu_t*d_m)/(pi*d_m - mu_t*p) + (F*mu_c*d_c/2) % Nm
T_L = (F*d_m/2)*(pi*mu_t*d_m - p)/(pi*d_m + mu_t*p) + (F*mu_c*d_c/2)

eff_R = F*p/2/pi/T_R
eff_L = F*p/2/pi/T_L

%% 6
clear
clc

tpi = 6;
p = 1/tpi; % in
mu = .15;
d_t = 3/4; % in
d_m = d_t - p/2;
d_c = 1;
T = 8*3.5; % lbfin

syms F
eqn = T == (F*d_m/2)*(p + pi*mu*d_m)/(pi*d_m - mu*p) + (F*mu*d_c/2); % lbfin
F = double(solve(eqn, F)) % lbf
