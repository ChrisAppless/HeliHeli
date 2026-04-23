clear; close all; clc;

sim_res = sim('helicopter_sim.slx');


%%
grndDist = sim_res.logsout.getElement('Ground Distance');
vel      = sim_res.logsout.getElement('Velocity');
refVel   = sim_res.logsout.getElement('Reference Velocity');
alt      = sim_res.logsout.getElement('Altitude');

hover_start = 100;

% Ground Distance
figure;
plot(grndDist.Values.Time, grndDist.Values.Data, 'b', 'LineWidth', 1.5);
xline(hover_start, 'k--', 'Hover Starts', 'LineWidth', 1)
xlabel('Time (s)');
ylabel('Ground Distance (m)');
title('Ground Distance');
grid on;
exportgraphics(gcf, 'plots/ground_distance.pdf', 'ContentType', 'vector');

% Velocity and Reference Velocity
figure; hold on;
plot(vel.Values.Time,    vel.Values.Data,    'b',   'LineWidth', 1.5);
plot(refVel.Values.Time, refVel.Values.Data, 'r--', 'LineWidth', 1.5);
hold off
xlabel('Time (s)');
ylabel('Velocity (kts)');
title('Velocity vs Reference Velocity');
legend('Velocity', 'Reference Velocity');
xline(hover_start, 'k--', 'Hover Starts', 'LineWidth', 1, 'HandleVisibility','off')
ylim([-5 45])
grid on;
exportgraphics(gcf, 'plots/velocity.pdf', 'ContentType', 'vector');

% Zoomed start in Velocity and Reference Velocity
figure; hold on;
plot(vel.Values.Time,    vel.Values.Data,    'b',   'LineWidth', 1.5);
plot(refVel.Values.Time, refVel.Values.Data, 'r--', 'LineWidth', 1.5);
hold off
xlabel('Time (s)');
ylabel('Velocity (kts)');
title('Velocity vs Reference Velocity');
legend('Velocity', 'Reference Velocity');
xline(hover_start, 'k--', 'Hover Starts', 'LineWidth', 1, 'HandleVisibility','off')
ylim([-5 45])
xlim([0 hover_start])
grid on;
exportgraphics(gcf, 'plots/velocity_zoomed_decel.pdf', 'ContentType', 'vector');

% Zoomed hov in Velocity and Reference Velocity
figure; hold on;
plot(vel.Values.Time,    vel.Values.Data,    'b',   'LineWidth', 1.5);
plot(refVel.Values.Time, refVel.Values.Data, 'r--', 'LineWidth', 1.5);
hold off
xlabel('Time (s)');
ylabel('Velocity (kts)');
title('Velocity vs Reference Velocity');
legend('Velocity', 'Reference Velocity');
xline(hover_start, 'k--', 'Hover Starts', 'LineWidth', 1, 'HandleVisibility','off')
xlim([hover_start 450])
grid on;
exportgraphics(gcf, 'plots/velocity_zoomed_hov.pdf', 'ContentType', 'vector');


% Altitude
figure;
plot(alt.Values.Time, alt.Values.Data, 'b', 'LineWidth', 1.5);
xline(hover_start, 'k--', 'Hover Starts', 'LineWidth', 1)
xlabel('Time (s)');
ylabel('Altitude (ft)');
title('Altitude');
grid on;
exportgraphics(gcf, 'plots/altitude.pdf', 'ContentType', 'vector');