function [heli, atm] = GetParams()

% Helicopter Parameters
heli.m = 4309;
heli.payload = 1257;
heli.R = 14.63/2;
heli.N = 2;
heli.RPS = 5.4;
heli.maxV = 56.67;
heli.c = 0.53;
heli.range = 511; %km
heli.hoverceiling = 1220;
heli.Ntail = 2;
heli.Rtail = 2.59/2;
heli.cgdist = 2.07;
heli.Iyy = 17259;
heli.FM = 0.75;
heli.k = 1.15;

heli.sigma = heli.N*heli.c/(pi*heli.R);
heli.W = heli.m*9.81;
heli.Omega = heli.RPS * 2*pi;
heli.A = pi * heli.R^2;
heli.Mu = heli.maxV / (heli.Omega * heli.R);

% Atmosphere Parameters
atm.rho = 1.225;
atm.T = 287.3;