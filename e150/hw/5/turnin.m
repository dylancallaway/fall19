dbstop if error
format long
close all
clear
clc

% Search bounds
k1 = 80e9;
u1 = 30e9;
o1 = 1e7;
K1 = 4.3;
k2_minus = k1;
k2_plus = 10*k1;
u2_minus = u1;
u2_plus = 10*u1;
K2_minus = K1;
K2_plus = 10*K1;
o2_minus = o1;
o2_plus = 10*o1;
v2_minus = 0;
v2_plus = 2/3;

f = @(k2, u2, o2, K2, v2) myMaterial(k2, u2, o2, K2, v2);

interval = [
    k2_minus, k2_plus;
    u2_minus, u2_plus;
    o2_minus, o2_plus;
    K2_minus, K2_plus;
    v2_minus, v2_plus
    ];

parents = 10;
TOL_GA = 0;
G = 5000;
S = 100;
dv = 5;

[PI1, Orig1, Lambda1, best_designs1] = myGenetic(f, interval, parents, TOL_GA, G, S, dv);
best_cost1 = PI1(G, 1)

[PI2, Orig2, Lambda2, best_designs2] = myGenetic_noparents(f, interval, parents, TOL_GA, G, S, dv);
best_cost2 = PI2(G, 1)

[PI3, Orig3, Lambda3, best_designs3] = myGenetic_mutations(f, interval, parents, TOL_GA, G, S, dv);
best_cost3 = PI3(G, 1)

% Deliverables
%% #2
figure(1);
loglog(PI1(:,1))
hold on
loglog(PI2(:,1))
loglog(PI3(:,1))
legend('Ten Parents', 'Zero Parents', 'Ten Parents, Child Mutations')
xlabel('Generation')
ylabel('Best Design Cost')
hold off

%% #3
figure(2);
subplot(1, 3, 1)
loglog(best_designs1(1, :)/10^9)
hold on
loglog(best_designs1(2, :)/10^9)
loglog(best_designs1(3, :)/10^6)
loglog(best_designs1(4, :))
loglog(best_designs1(5, :))
legend('\kappa_2 [GPa]', '\mu_2 [GPa]', '\sigma_{e, 2} [MS/m]', 'K_2 [W/m-K]', 'v_2 [m^3/m^3]')
xlabel('Generation')
ylabel('Design Parameter Value')
ylim([10^-1, 10^3])
hold off

subplot(1, 3, 2)
loglog(best_designs2(1, :)/10^9)
hold on
loglog(best_designs2(2, :)/10^9)
loglog(best_designs2(3, :)/10^6)
loglog(best_designs2(4, :))
loglog(best_designs2(5, :))
legend('\kappa_2 [GPa]', '\mu_2 [GPa]', '\sigma_{e, 2} [MS/m]', 'K_2 [W/m-K]', 'v_2 [m^3/m^3]')
xlabel('Generation')
ylabel('Design Parameter Value')
ylim([10^-1, 10^3])
hold off

subplot(1, 3, 3)
loglog(best_designs3(1, :)/10^9)
hold on
loglog(best_designs3(2, :)/10^9)
loglog(best_designs3(3, :)/10^6)
loglog(best_designs3(4, :))
loglog(best_designs3(5, :))
legend('\kappa_2 [GPa]', '\mu_2 [GPa]', '\sigma_{e, 2} [MS/m]', 'K_2 [W/m-K]', 'v_2 [m^3/m^3]')
xlabel('Generation')
ylabel('Design Parameter Value')
ylim([10^-1, 10^3])
hold off

%% #4
kdes = 111e9;
udes = 47e9;
odes = 2e7;
Kdes=6.2;

best_design(:, 1) = best_designs1(:, end);
best_design(:, 2) = best_designs2(:, end);
best_design(:, 3) = best_designs3(:, end);

k_diff(1) = (kdes - best_design(1, 1)) / kdes;
u_diff(1) = (udes - best_design(2, 1)) / udes;
o_diff(1) = (odes - best_design(3, 1)) / odes;
K_diff(1) = (Kdes - best_design(4, 1)) / Kdes;

k_diff(2) = (kdes - best_design(1, 2)) / kdes;
u_diff(2) = (udes - best_design(2, 2)) / udes;
o_diff(2) = (odes - best_design(3, 2)) / odes;
K_diff(2) = (Kdes - best_design(4, 2)) / Kdes;

k_diff(3) = (kdes - best_design(1, 3)) / kdes;
u_diff(3) = (udes - best_design(2, 3)) / udes;
o_diff(3) = (odes - best_design(3, 3)) / odes;
K_diff(3) = (Kdes - best_design(4, 3)) / Kdes;

type = ["Ten Parents", "Zero Parents", "Ten Parents, Child Mutations"];
param = ["\kappa_2", "\mu_2", "\sigma_{e, 2}", "K_2"];
for i = 1:3
    fprintf("%s & ", type(i));
    fprintf("%.3f & ", k_diff(i))
    fprintf("%.3f & ", u_diff(i))
    fprintf("%.3f & ", o_diff(i))
    fprintf("%.3f %s", K_diff(i), "\\")
    fprintf("\n")
end

%% #5
sum_diff(1) = sum([abs(k_diff(1)), abs(u_diff(1)), abs(o_diff(1)), abs(K_diff(1))]);
sum_diff(2) = sum([abs(k_diff(2)), abs(u_diff(2)), abs(o_diff(2)), abs(K_diff(2))]);
sum_diff(3) = sum([abs(k_diff(3)), abs(u_diff(3)), abs(o_diff(3)), abs(K_diff(3))]);

for i = 1:3
    fprintf("%s & ", type(i));
    fprintf("%.3f %s ", sum_diff(i), "\\")
    fprintf("\n")
end

%% #7
best_designs = {best_designs1, best_designs2, best_designs3};
for c = 1:3
    best_designs{c}(1, :) = best_designs{c}(1, :) / 10^9;
    best_designs{c}(2, :) = best_designs{c}(2, :) / 10^9;
    best_designs{c}(3, :) = best_designs{c}(3, :) / 10^6;
    for i = 1:4
        fprintf("%d & ", i);
        % nu
        best_designs{c}(6, end-i+1) = (3*best_designs{c}(1, end-i+1) - 2*best_designs{c}(2, end-i+1)) / (2*best_designs{c}(2, end-i+1) + 6*best_designs{c}(1, end-i+1));
        % E^y
        best_designs{c}(7, end-i+1) = 2*best_designs{c}(1, end-i+1)*(1+best_designs{c}(5, end-i+1));
        for j = 1:7
            if j == 7
                fprintf("%.3f %s", best_designs{c}(j, end-i+1), "\\")
            else
                fprintf("%.3f & ", best_designs{c}(j, end-i+1))
            end
        end
        fprintf("\n")
    end
    fprintf("\n")
end