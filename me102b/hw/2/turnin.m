%% HW2
%% Dylan Callaway

%% 1
close all
clear
clc

torque_req = 10*3/8;
torque_max = 29.5;
V_max = 15;
V_req = V_max*torque_req/torque_max % V

%% 2
clear
clc

K_T = 16.1;
R_coil = 1.33;
% Is R_coil here the total resistance of the coils, or the
% resistance of a single coil that we need to add in parallel to the
% other 2 (N)?
omega_NL = 10300;
V_op = 18; % Assume ideal voltage source/constant voltage.
torque_max = 24.2;
K_e = V_op/omega_NL;
% Assuming the question is asking what the rotational speed is when the
% max continuous torque is applied.
omega = omega_NL - R_coil/K_T/K_e*torque_max % rpm

%% 3
clear
clc

omega_NL = 11500;
torque_stall = 4.47;
torque = 1.5;
% Assume constant voltage.
omega_load = torque/torque_stall*omega_NL % RPM

%% 4
clear
clc

K_T = 1.657;
K_e = 1.23;
R = 20.3;
V = 14;
torque = .15;
dc0 = .25;
dc1 = .85;
V0 = dc0*V;
V1 = dc1*V;
omega0 = V0/K_e - R/K_e/K_T*torque % KRPM
omega1 = V1/K_e - R/K_e/K_T*torque % KRPM

%% 5
clear
clc

delay = 8;
T = 25;
Gp = 20/.1; % 20rpm / 10% dc
K_p = 1.2*T/delay/Gp
K_i = .5/delay
K_d = .5*delay