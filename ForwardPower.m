clear; close all; clc

Omega = 5.4;
R = 7.315;
V = linspace(0,56.67,100);

% Profile drag power
P_hover_BEM = 550.9e3;
Mu = V / (Omega * R);
P_profile = P_hover_BEM*(1+Mu.^2);

% Rotor drag power
