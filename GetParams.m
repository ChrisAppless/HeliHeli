function [heli, atm] = GetParams()

% Helicopter Parameters
heli.m = 4309;
heli.R = 14.63/2;
heli.N = 2;
heli.RPS = 5.4;
heli.maxV = 56.67;

heli.W = heli.m*9.81;
heli.omega = heli.RPS * 2*pi;
heli.A = pi * heli.R^2;

% Atmosphere Parameters
atm.rho = 1.225;