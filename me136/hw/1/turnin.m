%% 1
close all
clear
clc

m = 900;
g = 9.81;
rho = 1.225;
C_L = 1.6;
S = 16;
F_L = m*g;

v_min = sqrt( 2*F_L / rho / S / C_L ) % m/s

%% 2
clear
clc

% 1 A380
fprintf("\nA380\n")
v = 903*1000/60/60;
rho = .4135;
mu = 1.458E-5;
c = 10.6; % avg chord from http://www.dept.aoe.vt.edu/~mason/Mason_f/A380Dean.pdf
Re = v*c*rho/mu
v_sound = 295.2;
mach_number = v/v_sound
AR = 7.8
b = 79.6;
S = b^2/AR;
g = 9.776;
m_L = 650*1000;
f_L = m_L*g;
C_L = 2*f_L/rho/v^2/S
C_D = C_L^2/pi/AR

% 2 DG-808C
fprintf("\nDG-808C\n")
v = 150*1000/60/60;
rho = .7364;
mu = 1.628E-5;
c = .7; % avg chord from https://cevans.me/GLASS/Documentation/DG/DG800/images/Manual_800B.pdf
Re = v*c*rho/mu
v_sound = 320.5;
mach_number = v/v_sound
b = 7.06;
S = 11.8;
AR = b^2/S
g = 9.791;
m_L = 156;
f_L = m_L*g;
C_L = 2*f_L/rho/v^2/S
C_D = C_L^2/pi/AR

% 3 F-22
fprintf("\nF-22\n")
v = 1900*1000/60/60;
rho = .4135;
mu = 1.458E-5;
c = 10.5;
Re = v*c*rho/mu
v_sound = 295.2;
mach_number = v/v_sound
AR = 2.36
b = 13.56;
S = 78.04;
g = 9.776;
m_L = 38000;
f_L = m_L*g;
C_L = 2*f_L/rho/v^2/S
C_D = C_L^2/pi/AR

% 4 eBee
fprintf("\neBee\n")
v = 60*1000/60/60;
rho = 1.225;
mu = 1.789E-5;
c = .18; % approx chord from https://www.yanmu.co.id/2015/11/ebee-sensefly.html
Re = v*c*rho/mu
v_sound = 340.3;
mach_number = v/v_sound
AR = b^2/S
b = .96;
S = .5*.18*.96; % approx
g = 9.81;
m_L = .7; % approx
f_L = m_L*g;
C_L = 2*f_L/rho/v^2/S
C_D = C_L^2/pi/AR

%% 3
clear
clc

g = 9.81;
C_T = 6.41E-6;
gamma = .017;
m = .5;

% a
f_total = m*g;
f_req = f_total/4;

omega = sqrt( f_req/C_T ) % rad/s

% b
tau = gamma*f_req;
P = tau*omega % W

%% 4
clear
clc

% a
g = 9.81;
m = 1000;
b = 9;
S = .7*4.5 + 9*.9;
AR = b^2/S

% b
syms alpha
C_l = 2*pi*alpha;
C_L = AR*C_l/(AR+2)

% c
f_L = m*g;
rho = 1.225;
v = 55;
eqn = 2*f_L/rho/S/C_L == v;
alpha_req = double(solve(eqn, alpha)) % degrees

%% 5
clear
clc

g = 9.81;
m = 7;
d = .217;
rho = 1.225;
mu = 18.27E-6;
v = 0:150;
Re = v*d*rho/mu;
figure
plot(v, Re)
title("Re vs. v")
xlabel("v (m/s)")
ylabel("Re")

f_req = m*g;
S = pi*d^2/4;
C_D = .16;
% I plotted Re vs. v to see what range of Re we might be in, then solved
% with C_D=.4, C_D=.15 and C_D=.16 until the terminal velocity at a
% given C_D matched the C_D at the given Re.
v = sqrt( 2*f_req/rho/S/C_D )