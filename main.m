clear; close all; clc

PLOT_DIR = "plots/";
if ~exist(PLOT_DIR, 'dir'), mkdir(PLOT_DIR); end
set(0, 'DefaultFigureVisible', 'off');

[heli, atm] = GetParams;

perf = struct();

%% Induced velocity Hover
perf.vi_hov = sqrt(heli.W/(2*atm.rho*pi*heli.R^2));

%% Induced velocity forward flight
perf.Vbar_range = (0:1:heli.maxV)./perf.vi_hov;
perf.vi_ff = zeros(size(perf.Vbar_range));
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
    perf.vi_ff(i) = vi;
end

figure; hold on;
plot(perf.Vbar_range, perf.vi_ff, 'b-', 'LineWidth', 2);
plot(perf.Vbar_range, 1./perf.Vbar_range, 'k--')
plot((1./perf.vi_ff - perf.vi_ff), perf.vi_ff, 'k--'); hold off;
xlabel('Forward velocity V');
ylabel('Induced velocity v_i');
title('Induced velocity vs forward speed (normalized to v_{i,hov})');
grid on; ylim([0 1]); xlim([0 5]);
exportgraphics(gcf, PLOT_DIR + "vi_forward_flight.png");

%% Power calculation Hover
perf.P_hov_i = heli.W * perf.vi_hov;
perf.P_hov_ACT = perf.P_hov_i / heli.FM;

perf.P_hov_BEM =

%% Power calculation Forward flight

% Profile drag power
P_hover_BEM = 550.9e3;
perf.P_profile = P_hover_BEM*(1+heli.Mu.^2);

% Rotor drag power