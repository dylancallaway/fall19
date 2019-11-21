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

[PI1, Orig1, Lambda1, best_design1] = myGenetic(f, interval, parents, TOL_GA, G, S, dv);
best_cost1 = PI1(G, 1)

[PI2, Orig2, Lambda2, best_design2] = myGenetic_noparents(f, interval, parents, TOL_GA, G, S, dv);
best_cost2 = PI2(G, 1)

[PI3, Orig3, Lambda3, best_design3] = myGenetic_mutations(f, interval, parents, TOL_GA, G, S, dv);
best_cost3 = PI3(G, 1)

plot(PI1(:,1))
hold on
plot(PI2(:,1))
plot(PI3(:,1))
legend('normal', 'no parents', 'mutato')


%%
cost = f(best_design(1), best_design(2), best_design(3), best_design(4), best_design(5))


