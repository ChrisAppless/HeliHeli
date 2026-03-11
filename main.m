clear; close all; clc

[heli, atm] = GetParams;

%% Induced velocity Hover
vi_hov = sqrt(heli.W/(2*atm.rho*pi*heli.R^2));

%% Induced velocity forward flight
Vbar_range = (0:1:heli.maxV)./vi_hov;
vi_ff = zeros(size(Vbar_range));
alphad = 0;

for i = 1:length(Vbar_range)
    Vbar = Vbar_range(i);
    vi = 1; 
    for iter = 1:1000
        vi_new = 1 / sqrt((Vbar*cos(alphad))^2 + (Vbar*sin(alphad) + vi)^2);
        if abs(vi_new - vi) < 1e-6
            break;
        end
        vi = vi_new;
    end
    vi_ff(i) = vi;
end

figure; hold on;
plot(Vbar_range, vi_ff, 'b-', 'LineWidth', 2);
plot(Vbar_range, 1./Vbar_range, 'k--')
plot((1./vi_ff - vi_ff), vi_ff, 'k--'); hold off;
xlabel('Forward velocity V');
ylabel('Induced velocity v_i');
title('Induced velocity vs forward speed (normalized to v_{i,hov})');
grid on;
ylim([0 1]);
xlim([0 5]);

%%Forward flight power calculations

% Profile drag power
P_hover_BEM = 550.9e3;
Mu = heli.maxV / (heli.Omega * heli.R);
P_profile = P_hover_BEM*(1+Mu.^2);

% Rotor drag power
