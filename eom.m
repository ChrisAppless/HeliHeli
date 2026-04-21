function xdot = eom(x, theta_c, theta_0)
xdot = zeros(4,1);
% x = [u; w; q; tf] — W6_L11_P1 slide 9

% Parameters
m       = 4309;       R     = 7.315;      N       = 2;
c       = 0.53;       Omega = 33.9292;    W       = 42271.29;
Iyy     = 17259;      cgdist= 2.07;       Afront  = 4.46;
fuse_CD = 0.2;        m_blade = 101;      rho     = 1.225;

a_lift = 2*pi;
sigma  = N * c / (pi * R);
I_bl   = (1/3) * m_blade * R^2;
gamma  = rho * a_lift * c * R^4 / I_bl;

% Step 1: kinematics
V       = max(sqrt(x(1)^2 + x(2)^2), 0.01);
alpha_c = theta_c - atan2(x(2), x(1));
if x(1) < 0, alpha_c = alpha_c + pi; end
mu       = (V / (Omega * R)) * cos(alpha_c);
lambda_c = (V / (Omega * R)) * sin(alpha_c);
VOR      = V / (Omega * R);

% Step 2: iterate for lambda_i
li = sqrt(W / (2 * rho * pi * R^2)) / (Omega * R);
for iter = 1:100
    a1_temp    = ((8/3)*mu*theta_0 - 2*mu*(lambda_c+li) - (16/gamma)*(x(3)/Omega)) / (1 - 0.5*mu^2);
    CT_BEM     = (a_lift*sigma/4) * ((2/3)*theta_0*(1+1.5*mu^2) - (li+lambda_c));
    CT_Glauert = 2*li * sqrt((VOR*cos(alpha_c-a1_temp))^2 + (VOR*sin(alpha_c-a1_temp)+li)^2);
    li_new     = CT_BEM / (2 * sqrt((VOR*cos(alpha_c-a1_temp))^2 + (VOR*sin(alpha_c-a1_temp)+li)^2));
    if abs(li_new - li) < 1e-8, break; end
    li = li_new;
end

% Step 3: final a1, T, D
a1 = ((8/3)*mu*theta_0 - 2*mu*(lambda_c+li) - (16/gamma)*(x(3)/Omega)) / (1 - 0.5*mu^2);
CT = (a_lift*sigma/4) * ((2/3)*theta_0*(1+1.5*mu^2) - (li+lambda_c));
T  = CT * rho * (Omega*R)^2 * pi * R^2;
D  = Afront * fuse_CD * 0.5 * rho * V^2;

% Step 4: EOM
disc_tilt = theta_c - a1;
xdot = [-9.81*sin(x(4)) - (D/m)*(x(1)/V) + (T/m)*sin(disc_tilt) - x(3)*x(2);
         9.81*cos(x(4)) - (D/m)*(x(2)/V) - (T/m)*cos(disc_tilt) + x(3)*x(1);
        -(T*cgdist/Iyy)*sin(disc_tilt);
         x(3)];
end
