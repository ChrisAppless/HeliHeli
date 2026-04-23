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

% Question 1
theta_0 = deg2rad(6);   % collective pitch [rad]
A1      = deg2rad(2);   % longitudinal cyclic [rad]
B1      = deg2rad(1);   % lateral cyclic [rad]

V = heli.forV;                  % [m/s]
q = deg2rad(heli.q);    % pitch rate [rad/s]
p = deg2rad(heli.p);    % roll rate [rad/s]

% Lock number
a_lift = 2*pi;
I_bl   = (1/3) * heli.m_blade * heli.R^2;
gamma  = atm.rho * a_lift * heli.c * heli.R^4 / I_bl;

% Control plane angle of attack (level flight trim: rotor tilts forward to overcome drag)
D_para   = heli.Afront * heli.fuse_CD * 0.5 * atm.rho * V^2;
alpha_c  = -atan(D_para / heli.W);  % negative = nose down

% Non-dimensional inflow
mu       = V * cos(alpha_c) / (heli.Omega * heli.R);
lambda_c = V * sin(alpha_c) / (heli.Omega * heli.R);
lambda_i = vi_interp(V) / (heli.Omega * heli.R);

% Flapping coefficients (order matters: a0 first)
a0 = gamma/8 * (theta_0*(1 + mu^2) - (4/3)*(lambda_i + lambda_c));

a1 = (A1 - (16/gamma)*(q/heli.Omega) + (8/3)*mu*theta_0 - 2*mu*(lambda_c + lambda_i)) ...
     / (1 - mu^2/2);

b1 = (B1 - p/heli.Omega + (4/3)*mu*a0) ...
     / (1 + mu^2/2);

psi          = linspace(0, 360, 360);
flapping_ang = rad2deg(a0 - a1*cosd(psi) - b1*sind(psi));

figure;
plot(psi, flapping_ang, 'b-', 'LineWidth', 2);
xlabel('\psi (deg)');
ylabel('\beta (deg)');
title('Blade flapping angle over one revolution');
xlim([0 360])
grid on;
exportgraphics(gcf, PLOT_DIR + "P2_flapping_angle.png");

% Question 2 
xy_lin = linspace(-1, 1, 400);
[Xq, Yq] = meshgrid(xy_lin, xy_lin);

r_cart   = sqrt(Xq.^2 + Yq.^2) * heli.R;
psi_cart = atan2d(Xq, -Yq);
psi_cart(psi_cart < 0) = psi_cart(psi_cart < 0) + 360;

theta_cart    = theta_0 - A1*cosd(psi_cart) - B1*sind(psi_cart);
beta_cart     = a0 - a1*cosd(psi_cart) - b1*sind(psi_cart);
beta_dot_cart = heli.Omega * (a1*sind(psi_cart) - b1*cosd(psi_cart));

V_perp_cart = V*sin(alpha_c) + vi_interp(V) ...
            + beta_dot_cart.*r_cart ...
            - q*r_cart.*cosd(psi_cart) ...
            + V*cos(alpha_c)*cosd(psi_cart).*beta_cart;

V_tan_cart = heli.Omega*r_cart + V*cos(alpha_c)*sind(psi_cart);

alpha_cart = rad2deg(theta_cart - V_perp_cart ./ V_tan_cart);
alpha_cart(Xq.^2 + Yq.^2 > 1) = NaN;             % outside disc
alpha_cart(Xq.^2 + Yq.^2 + 2*mu*Xq <= 0) = NaN;  % reverse flow region

figure;
contour(Xq, Yq, alpha_cart, -2:0.5:7, 'ShowText', 'on', 'LineWidth', 1.2);
xlabel('Blade position (nondimensional)');
ylabel('Blade position (nondimensional)');
title('Curves of constant angle of attack');
text(-0.7, -1.15, 'Retreating side', 'HorizontalAlignment', 'center');
text( 0.7, -1.15, 'Advancing side',  'HorizontalAlignment', 'center');
axis equal; grid on; colorbar;
exportgraphics(gcf, PLOT_DIR + "P2_AoA_disc.png");

% Question 3
fprintf("$a_0$ & %.4f° \\\\\n", rad2deg(a0))
fprintf("$a_1$ & %.4f° \\\\\n", rad2deg(a1))
fprintf("$b_1$ & %.4f° \\\\\n", rad2deg(b1))



%% Part 3 - Phugoid Analysis

% Non-dimensionalization factors from Pavel's nomenclature
A_force = atm.rho * pi*heli.R^2 * (heli.Omega*heli.R);          % for X_u, Z_u (δF/δvel)
A_moment = atm.rho * pi*heli.R^2 * heli.Omega * heli.R^3;        % for M_q (δM/δangvel)
A_moment_vel = atm.rho * pi*heli.R^2 * heli.R * heli.Omega;     % for M_u (δM/δvel)

% Hover trim values needed for derivatives
C_T   = heli.W / (atm.rho * pi*heli.R^2 * (heli.Omega*heli.R)^2);
lam_i0_hov = sqrt(C_T/2);
lam_i_hov  = (1 + heli.k_inflow) * lam_i0_hov;  % heli.k_inflow = 0.15 typically
lam_hov    = -lam_i_hov;                         % alpha = 0 at hover
theta0_hov = (6*C_T/(heli.sigma*a_lift) - 1.5*lam_hov);  % from closed-form hover collective

% Non-dimensional CG offsets
h_bar = heli.cgdist / heli.R;      % vertical CG-to-hub distance (nondim)
f_bar = 0 / heli.R;      % horizontal CG-to-hub offset (nondim)

% Hub moment coefficient (flap hinge / spring contribution)
% For articulated: C_mE = e_bar * N_b * M_beta / (2 * rho * A * (Omega*R)^2 * R)
% For hingeless: uses equivalent flap spring K_beta
A     = pi * heli.R^2;
C_mE = heli.e_bar * heli.m_blade / (atm.rho * A * heli.R);

% Flapping derivative at hover
da1_dmu_hov = (8/3)*theta0_hov + 2*lam_hov;

% H-force derivative at hover  
C_dm = getCD(perf.hov.alpha);  % mean drag coefficient
lam_D_hov = lam_hov;            % at hover these coincide
dCHD_dmu_hov = heli.sigma*C_dm/4 - (a_lift*heli.sigma*lam_D_hov/8)*da1_dmu_hov;

% Flap response to pitch rate at hover
da1_dqbar_hov = -16/gamma;

% Non-dimensional stability derivatives (Pavel's lowercase)
x_u = -C_T * da1_dmu_hov - dCHD_dmu_hov;
m_u_prime = C_T*h_bar*da1_dmu_hov + C_mE*da1_dmu_hov;  % f_bar * dCT/dmu term vanishes at hover
m_q_prime = -(16/gamma) * (C_T*h_bar + C_mE);

% Convert to dimensional form for the characteristic equation
I_y = heli.Iyy;
m   = heli.m;

X_u = x_u * A_force;
M_u = m_u_prime * A_moment_vel * heli.R;   
M_q = m_q_prime * A_moment;

% Dimensional derivatives per unit mass / inertia
Xu_hat = X_u / m;
Mu_hat = M_u / I_y;
Mq_hat = M_q / I_y;

% Build the state matrix for [u; theta; q]
A_phugoid = [Xu_hat,  -9.81,  0;
             0,        0,     1;
             Mu_hat,   0,     Mq_hat];

% Find eigenvalues (roots of characteristic polynomial)
eigs_phugoid = eig(A_phugoid);

% Or equivalently via polynomial roots
char_poly = [1, -(Xu_hat + Mq_hat), Xu_hat*Mq_hat, 9.81*Mu_hat];
roots_poly = roots(char_poly);

% Report
fprintf("\n--- Phugoid Analysis (Hover) ---\n");
fprintf("X_u = %.4e  [1/s]\n", Xu_hat);
fprintf("M_u = %.4e  [rad/(m·s)]\n", Mu_hat);
fprintf("M_q = %.4e  [1/s]\n", Mq_hat);
fprintf("\nEigenvalues:\n");
for k = 1:length(eigs_phugoid)
    lam = eigs_phugoid(k);
    if abs(imag(lam)) > 1e-6
        omega_n = abs(lam);
        zeta    = -real(lam)/omega_n;
        period  = 2*pi/abs(imag(lam));
        fprintf("  %+.4f %+.4fi   omega_n=%.3f rad/s, zeta=%.3f, T=%.2f s\n", ...
                real(lam), imag(lam), omega_n, zeta, period);
    else
        t_half = log(2)/abs(real(lam));
        fprintf("  %+.4f            real mode, t_{1/2 or 2}=%.2f s\n", real(lam), t_half);
    end
end

% Plot poles on complex plane
figure; hold on;
plot(real(eigs_phugoid), imag(eigs_phugoid), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xline(0, 'k--'); yline(0, 'k--');
xlabel('Real part (1/s)');
ylabel('Imaginary part (1/s)');
title('Phugoid poles in hover');
grid on; axis equal;
exportgraphics(gcf, PLOT_DIR + "P3_phugoid_poles.png");