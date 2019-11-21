function [cost] = myMaterial(k2, u2, o2, K2, v2)
k1 = 80e9;
u1 = 30e9;
o1 = 1e7;
K1 = 4.3;
kdes = 111e9;
udes = 47e9;
odes = 2e7;
Kdes=6.2;
TOLk = 0.5;
TOLu = 0.5;
TOLK = 0.5;
TOLo = 0.8;
W1 = 1/3;
W2 = 1/3;
W3 = 1/3;
w1 = 1;
w2 = 2;
phi = 0.5;

v1 = 1 - v2;

% Bounds and effective properties
% Bulk modulus
kminus =  k1 + v2 / ((1/(k2 - k1)) + (3*(1-v2)/(3*k1 + 4*u1)));
kplus = k2 + (1 - v2) / ((1/(k1 - k2)) + 3*v2/(3*k2 + 4*u2));
k = phi*kminus + (1-phi)*kplus;

% Shear modulus
uminus = u1 + v2 / (1/(u2 - u1) +...
    6*(1 - v2)*(k1 + 2*u1)/5/u1/(3*k1 + 4*u1));
uplus = u2 + (1 - v2) / (1/(u1 - u2) +...
    6*v2*(k2 + 2*u2)/5/u2/(3*k2 + 4*u2));
u = phi*uminus + (1-phi)*uplus;

% Electrical conductivity
ominus = o1 + v2 / (1/(o2 - o1) + (1 - v2)/3/o1);
oplus = o2 + (1 - v2)/(1/(o1-o2) + v2/3/o2);
o = phi*ominus + (1-phi)*oplus;

% Thermal conductivity
Kminus = K1 + v2 / (1/(K2 - K1) + (1 - v2)/3/K1);
Kplus = K2 + (1 - v2) / (1/(K1 - K2) + v2/3/K2);
K = phi*Kminus + (1-phi)*Kplus;

% Joule heating
CE1 = (o2 - o) / (o2-o1) / (1 - v2);
CE2 = (o - o1) / (o2-o1) / v2;

CJ1 = o1 / (o*(1-v2)) * (o2 - o) / (o2 - o1);
CJ2 = o2 / o / v2 * (o - o1) / (o2 - o1);

% Thermal load shares
Ctheta2 = (K - K1)/v2/(K2 - K1);
Cq2 = K2*Ctheta2/K;
Cq1 = (1 - v2*Cq2)/(1 - v2);

% Local fields
% Stress
Co2k = k2*(k - k1)/v2/k/(k2 - k1);
Co2u = u2*(u - u1)/v2/u/(u2 - u1);
Co1k = (1 - v2*Co2k)/v1;
Co1u = (1 - v2*Co2u)/v1;

% Cost functions
% Electrical
w2h = ~(CJ1*CE1-TOLo)/TOLo <= 0;
w3h = ~(CJ2*CE2-TOLo) <= 0;
cost_elec = w1*abs((odes-o)/odes) + w2h*abs((CJ1*CE1-TOLo)/TOLo)...
    + w3h*abs((CJ2*CE2-TOLo)/TOLo);

% Thermal
w2h = ~(Cq1-TOLK)/TOLK <= 0;
w3h = ~(Cq2-TOLK)/TOLK <= 0;
cost_thermal = w1*abs((Kdes-K)/Kdes) + w2h*abs((Cq1-TOLK)/TOLK)...
    + w3h*abs((Cq2-TOLK)/TOLK);

% Mechanical
w3h = ~Co2k <= TOLk;
w4h = ~Co2u <= TOLu;
w5h = ~Co1k <= TOLk;
w6h = ~Co1u <= TOLu;
cost_mech = w1*abs((kdes-k)/kdes) + w2*abs((udes-u)/udes)...
    + w3h*abs((Co2k-TOLk)/TOLk) + w4h*abs((Co2k-TOLu)/TOLu)...
    + w5h*abs((Co1k-TOLk)/TOLk) + w6h*abs((Co1u-TOLu)/TOLu);

% Total
cost = W1*cost_elec + W2*cost_thermal + W3*cost_mech;
end

