clear; close all; clc

PLOT_DIR = "plots/";
if ~exist(PLOT_DIR, 'dir'), mkdir(PLOT_DIR); end
%set(0, 'DefaultFigureVisible', 'off');

naca0012_polar = readmatrix("NACA0012Polar.txt");
getCD = @(alpha) interp1(naca0012_polar(:,1), naca0012_polar(:,4), alpha, 'linear');

[heli, atm] = GetParams;

perf = struct();

%% Part 1

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
perf.hov.alpha = 6*heli.W ./ (heli.N .* atm.rho .* heli.Omega.^2 .* heli.R.^3 .* heli.c .* 2*pi *(pi/180));

perf.hov.P_ideal = heli.W * perf.hov.vi;
perf.hov.P_ACT = perf.hov.P_ideal / heli.FM;

perf.hov.P_i  = heli.k * heli.W * perf.hov.vi;
perf.hov.P_p  = (heli.sigma * getCD(perf.hov.alpha) / 8) * atm.rho * (heli.Omega * heli.R)^3 * pi * heli.R^2;
perf.hov.P_BEM = perf.hov.P_i + perf.hov.P_p;

fprintf("Ideal Power required to hover: %.3fkW\n", perf.hov.P_ideal.*1e-3)
fprintf("ACT Power required to hover: %.3fkW\n", perf.hov.P_ACT.*1e-3)
fprintf("BEM Power required to hover: %.3fkW\n", perf.hov.P_BEM.*1e-3)

%% Power calculation Forward flight
vi_interp = @(V) perf.hov.vi .* interp1(perf.Vbar_range, perf.ff.vi, V ./ perf.hov.vi, "linear", "extrap");
getMu = @(V) V ./ (heli.Omega .* heli.R);

% Profile drag + rotor drag power (Bennett approximation, eq 14)
getPp = @(V) perf.hov.P_p*(1+4.65.*getMu(V).^2);

% Induced power 
getIP = @(V) heli.k*heli.m*9.81*vi_interp(V);

% Parasite drag power
getPpara = @(V) (heli.Afront .* heli.fuse_CD) .* 0.5 * atm.rho * V.^3;

% Tail rotor induced drag
T_tr = @(V) (getPp(V) + getIP(V)) ./ heli.Omega ./ heli.tailarm;
perf.hov.tr.vi = sqrt(T_tr(0)./(2.*atm.rho*pi*heli.Rtail^2));
perf.hov.tr.alpha = 6*T_tr(0) ./ (heli.Ntail .* atm.rho .* heli.OmegaTail.^2 .* heli.Rtail.^3 .* heli.ctail .* 2*pi *(pi/180));
P_i_tr = @(V) 1.1.*heli.k_tr .* T_tr(V) .* perf.hov.tr.vi;
P_p_tr = @(V) (heli.sigma_tr .* getCD(perf.hov.tr.alpha) ./ 8) * atm.rho * (heli.OmegaTail * heli.Rtail)^3 * pi * heli.Rtail^2;
getPBEMtr = @(V) P_i_tr(V) + P_p_tr(V);

% Total Power
P_tot = @(V) getPp(V) + getIP(V) + getPpara(V) + getPBEMtr(V);

V_ext = heli.maxV + 40;
figure; hold on;
% Solid lines within Vmax
fplot(@(V) P_tot(V)*1e-3,    [0 heli.maxV], 'k-',  'LineWidth', 2);
fplot(@(V) getPp(V)*1e-3,    [0 heli.maxV], 'b-');
fplot(@(V) getIP(V)*1e-3,    [0 heli.maxV], 'r-');
fplot(@(V) getPpara(V)*1e-3, [0 heli.maxV], 'g-');
fplot(@(V) getPBEMtr(V)*1e-3,[0 heli.maxV], 'm-');
% Dashed lines beyond Vmax
fplot(@(V) P_tot(V)*1e-3,    [heli.maxV V_ext], 'k--', 'LineWidth', 2);
fplot(@(V) getPp(V)*1e-3,    [heli.maxV V_ext], 'b--');
fplot(@(V) getIP(V)*1e-3,    [heli.maxV V_ext], 'r--');
fplot(@(V) getPpara(V)*1e-3, [heli.maxV V_ext], 'g--');
fplot(@(V) getPBEMtr(V)*1e-3,[heli.maxV V_ext], 'm--');
% Vertical line at Vmax
xline(heli.maxV, 'k:', 'LineWidth', 1.5, 'Label', sprintf('V_{max}=%.1f m/s', heli.maxV));
% Minimum power point
[V_opt, P_min] = fminbnd(@(V) P_tot(V), 1, heli.maxV);
plot(V_opt, P_min*1e-3, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(V_opt, P_min*1e-3+20, sprintf('V=%.1f m/s\nP=%.1f kW', V_opt, P_min*1e-3), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
hold off;
legend('P_{total}', 'P_{profile}', 'P_{induced}', 'P_{parasite}', 'P_{tail rotor}', '', '', '', '', '', '', 'P_{min}', 'Location', 'northeast');
xlabel('Forward velocity V (m/s)');
ylabel('Power (kW)');
title('Power vs forward velocity');
grid on;
exportgraphics(gcf, PLOT_DIR + "P_forward_flight.png");

%% Part 2

