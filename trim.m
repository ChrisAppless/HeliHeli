function [theta_0, theta_c] = trim(V)

[heli, atm] = GetParams;

a_lift = 2*pi;
sigma  = heli.N * heli.c / (pi * heli.R);

D_fus   = heli.Afront * heli.fuse_CD * 0.5 * atm.rho * V^2;
T       = sqrt(heli.W^2 + D_fus^2);                                          % W6_L11_P3 slide 3, eq 1
CT      = T / (atm.rho * (heli.Omega * heli.R)^2 * pi * heli.R^2);          % W6_L11_P2 slide 7, eq 12
alpha_d = atan(D_fus / heli.W);                                               % W6_L11_P3 slide 6, eq 8
mu      = V / (heli.Omega * heli.R);

glauert  = @(li) 2*li * sqrt((mu*sin(alpha_d) + li)^2 + (mu*cos(alpha_d))^2) - CT;  % W6_L11_P2 slide 10, eq 20
lambda_i = fzero(glauert, sqrt(CT/2), optimset('Display','off'));

x0  = [deg2rad(8); deg2rad(2)];
sol = fsolve(@(x) residuals(x, CT, lambda_i, mu, alpha_d, sigma, a_lift), x0, ...
             optimset('Display','off'));

theta_0 = sol(1);
theta_c = sol(2);

end

function res = residuals(x, CT, lambda_i, mu, alpha_d, sigma, a_lift)
    theta_0  = x(1);
    theta_c  = x(2);
    lambda_c = mu * (alpha_d + theta_c);  % W6_L11_P3 slide 8, eq 12

    % f: CT_elem = CT  (W6_L11_P2 slide 8, eq 16 / W6_L11_P3 slide 8, eq 11)
    CT_elem = (a_lift * sigma / 4) * ((2/3)*theta_0*(1 + 1.5*mu^2) - (lambda_i + lambda_c));
    f = CT_elem - CT;

    % g: theta_c = a1 in trim  (W6_L11_P3 slide 9)
    a1_flap = ((8/3)*mu*theta_0 - 2*mu*(lambda_i + lambda_c)) / (1 - 0.5*mu^2);
    g = theta_c - a1_flap;

    res = [f; g];
end