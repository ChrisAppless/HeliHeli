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

V2       = heli.forV;
q_rate   = deg2rad(heli.q);
p_rate   = deg2rad(heli.p);
theta_0  = deg2rad(6);
theta_1c = deg2rad(1);
theta_1s = deg2rad(2);

heli.m_blade = 100;
heli.Ib      = (1/3) * heli.m_blade * heli.R^2;
Cla          = 2*pi;
gamma_lock   = atm.rho * Cla * heli.c * heli.R^4 / heli.Ib;

alphad = 0;
mu     = V2 * cos(alphad) / (heli.Omega * heli.R);
vi_V2  = vi_interp(V2);
lambda = (V2 * sin(alphad) + vi_V2) / (heli.Omega * heli.R);

q_bar = q_rate / heli.Omega;
p_bar = p_rate / heli.Omega;

theta_fn = @(psi) theta_0 + theta_1c*cos(psi) + theta_1s*sin(psi);
uT_fn    = @(x, psi) x + mu*sin(psi);
uP_fn    = @(x, psi, beta, bdot) lambda + x.*bdot + mu*beta.*cos(psi) ...
           + x.*(p_bar*sin(psi) - q_bar*cos(psi));

g = gamma_lock;
M = [ 1,         0,              0;
      0,         g*(2-mu^2)/16,  0;
     -g*mu/6,    0,              g*(2+mu^2)/16 ];
r = [ g*theta_0*(1+mu^2)/8 + g*mu*theta_1s/6 - g*lambda/6 - g*mu*p_bar/12;
      g*mu*theta_0/3 + g*theta_1s*(2+3*mu^2)/16 - g*lambda*mu/4 - g*p_bar/8 + 2*q_bar;
     -g*theta_1c*(2+mu^2)/16 - g*q_bar/8 - 2*p_bar ];
coef = M \ r;
a0 = coef(1);
a1 = coef(2);
b1 = coef(3);

psi_plot  = linspace(0, 2*pi, 361)';
beta_plot = a0 - a1*cos(psi_plot) - b1*sin(psi_plot);
bdot_plot = a1*sin(psi_plot) - b1*cos(psi_plot);

perf.p2.a0 = rad2deg(a0);
perf.p2.a1 = rad2deg(a1);
perf.p2.b1 = rad2deg(b1);
perf.p2.gamma  = gamma_lock;
perf.p2.mu     = mu;
perf.p2.lambda = lambda;

fprintf('\n=== Part 2 - Rotor Dynamics ===\n');
fprintf('Lock number gamma        = %.3f\n', gamma_lock);
fprintf('Advance ratio mu         = %.4f\n', mu);
fprintf('Inflow ratio lambda      = %.4f\n', lambda);
fprintf('Coning a0                = %.3f deg\n', rad2deg(a0));
fprintf('Longitudinal tilt a1     = %.3f deg\n', rad2deg(a1));
fprintf('Lateral tilt      b1     = %.3f deg\n', rad2deg(b1));

figure;
plot(rad2deg(psi_plot), rad2deg(beta_plot), 'b-', 'LineWidth', 2);
xline(90,  'k:', 'advancing');
xline(180, 'k:', 'front');
xline(270, 'k:', 'retreating');
xlabel('\psi (deg)'); ylabel('\beta (deg)');
title(sprintf('Blade flapping  V=%g m/s, q=%g°/s, p=%g°/s', V2, heli.q, heli.p));
grid on; xlim([0 360]);
exportgraphics(gcf, PLOT_DIR + "P2_flapping_angle.png");

x_075     = 0.75;
uT_075    = uT_fn(x_075, psi_plot);
uP_075    = uP_fn(x_075, psi_plot, beta_plot, bdot_plot);
alpha_075 = theta_fn(psi_plot) - atan2(uP_075, uT_075);

figure;
plot(rad2deg(psi_plot), rad2deg(alpha_075), 'LineWidth', 2);
xline(90,  'k:', 'advancing');
xline(270, 'k:', 'retreating');
xlabel('\psi (deg)'); ylabel('\alpha (deg)');
title('Angle of attack at r/R = 0.75');
grid on; xlim([0 360]);
exportgraphics(gcf, PLOT_DIR + "P2_AoA_psi.png");

x_grid    = linspace(0.15, 1, 60);
[PSI, X]  = meshgrid(psi_plot, x_grid);
BETA      = repmat(beta_plot', length(x_grid), 1);
BDOT      = repmat(bdot_plot', length(x_grid), 1);
THETA     = theta_0 + theta_1c*cos(PSI) + theta_1s*sin(PSI);
UT_g      = X + mu*sin(PSI);
UP_g      = lambda + X.*BDOT + mu*BETA.*cos(PSI) + X.*(p_bar*sin(PSI) - q_bar*cos(PSI));
ALPHA_deg = rad2deg(THETA - atan2(UP_g, UT_g));

Xcart =  X .* sin(PSI);
Ycart = -X .* cos(PSI);

figure; hold on;
contourf(Xcart, Ycart, ALPHA_deg, 20, 'LineColor', 'none');
contour (Xcart, Ycart, ALPHA_deg, 12, 'k', 'LineWidth', 0.4, 'ShowText', 'on');
plot(cos(linspace(0,2*pi,200)), sin(linspace(0,2*pi,200)), 'k-', 'LineWidth', 1.5);
axis equal; colorbar; grid on;
xlabel('\leftarrow Retreating    Advancing \rightarrow');
ylabel('\leftarrow Aft    Forward \rightarrow');
title('\alpha on rotor disc [deg]');
xlim([-1.1 1.1]); ylim([-1.1 1.1]); hold off;
exportgraphics(gcf, PLOT_DIR + "P2_AoA_disc.png");

