function [heli, atm] = GetParams()

% Helicopter Parameters
heli.m = 4309;
heli.payload = 1257;
heli.R = 14.63/2;
heli.N = 2;
heli.RPS = 5.4;
heli.maxV = 56.67;
heli.forV = 20;
heli.q = 20; %Body pitch rate deg/s
heli.p = 10; %Roll rate deg/s
heli.c = 0.53;
heli.ctail = 0.213;
heli.range = 511; %km
heli.hoverceiling = 1220;
heli.Ntail = 2;
heli.Rtail = 2.59/2;
heli.cgdist = 2.07;
heli.tailarm = 8.53; % Horizontal distance of tail to cg
heli.Iyy = 17259;
heli.FM = 0.75;
heli.k = 1.15;
heli.CDp = 0.00322; % from airfoilTools 
heli.Afront = 4.46;
heli.fuse_CD = 0.2;
heli.k_tr = 1.4;
heli.OmegaTail = 174.13;
heli.m_blade = 101;     % blade mass [kg]
heli.k_inflow = 0.15;
heli.e_bar = 0.04;

heli.sigma = heli.Ntail*heli.ctail/(pi*heli.Rtail);
heli.sigma_tr = heli.N*heli.c/(pi*heli.R);
heli.W = heli.m*9.81;
heli.Omega = heli.RPS * 2*pi;
heli.A = pi * heli.R^2;
heli.Mu = heli.maxV / (heli.Omega * heli.R);

% Atmosphere Parameters
atm.rho = 1.225;
atm.T = 287.3;