clear; close all; clc

PLOT_DIR = "plots/";
if ~exist(PLOT_DIR, 'dir'), mkdir(PLOT_DIR); end
set(0, 'DefaultFigureVisible', 'off');

[heli, atm] = GetParams;

perf = struct();

%% Induced velocity Hover
perf.hov.vi = sqrt(heli.W/(2*atm.rho*pi*heli.R^2));

%% Induced velocity forward flight
perf.Vbar_range = (0:1:heli.maxV)./perf.hov.vi;
perf.ff.vi = zeros(size(perf.Vbar_range));
alphad = 0;

for i = 1:length(perf.Vbar_range)
    Vbar = perf.Vbar_range(i);
    vi = 1;
    for iter = 1:1000
        vi_new = 1 / sqrt((Vbar*cos(alphad))^2 + (Vbar*sin(alphad) + vi)^2);
        if abs(vi_new - vi) < 1e-6
            break;
        end
        vi = vi_new;
    end
    perf.ff.vi(i) = vi;
end

figure; hold on;
plot(perf.Vbar_range, perf.ff.vi, 'b-', 'LineWidth', 2);
plot(perf.Vbar_range, 1./perf.Vbar_range, 'k--')
plot((1./perf.ff.vi - perf.ff.vi), perf.ff.vi, 'k--'); hold off;
xlabel('Forward velocity V');
ylabel('Induced velocity v_i');
title('Induced velocity vs forward speed (normalized to v_{i,hov})');
grid on; ylim([0 1]); xlim([0 5]);
exportgraphics(gcf, PLOT_DIR + "vi_forward_flight.png");

%% Power calculation Hover
perf.hov.P_ideal = heli.W * perf.hov.vi;
perf.hov.P_ACT = perf.hov.P_ideal / heli.FM;

perf.hov.P_i  = heli.k * heli.W * perf.hov.vi;
perf.hov.P_p  = (heli.sigma * heli.CDp / 8) * atm.rho * (heli.Omega * heli.R)^3 * pi * heli.R^2;
perf.hov.P_BEM = perf.hov.P_i + perf.hov.P_p;

fprintf("Ideal Power required to hover: %.3fkW\n", perf.hov.P_ideal.*1e-3)
fprintf("ACT Power required to hover: %.3fkW\n", perf.hov.P_ACT.*1e-3)
fprintf("BEM Power required to hover: %.3fkW\n", perf.hov.P_BEM.*1e-3)

%% Power calculation Forward flight

% Profile drag power
perf.P_profile = perf.hov.P_BEM*(1+heli.Mu.^2);

% Rotor drag power
perf.P_drag = heli.sigma*heli.CDp/4 * atm.rho*(heli.Omega*heli.R)^3*pi*heli.R*heli.sigma^2;

% Induced power 
perf.P_induced = heli.k*heli.m*9.81*vi;

% Parasite drag power
perf.P_parasite = 0.5 * atm.rho * heli.maxV^3;