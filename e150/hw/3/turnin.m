dbstop if error
format compact
close all
clear
clc

%% Part 1
if 0
    theta_dot = [.2, -.2, 10];
    delta_v_d = [0, -1.2, 0];
    cost = myPrinter1(theta_dot, delta_v_d)
end

%% Part 2
if 1
    theta_dot = [.2, -.2, 10];
    delta_v_d = [0, -1.2, 0];
    cost = myPrinter2(theta_dot, delta_v_d)
end