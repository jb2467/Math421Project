load("Data/T.mat", "T");
load("Data/P.mat", "P");
load("Data/Rho.mat", "Rho");
load("Data/TemperatureL96.mat");
load("Data/particle_velocityL96.mat");
load("Data/U.mat", "U");

Nx = 100;
Nt = 2000;

Temperature = reshape(Temperature, [Nx, 1, 1, Nt]);
particle_velocity = reshape(particle_velocity, [Nx, 1, 1, Nt]);

TemperatureDim = 3*Temperature + T;
Ux = U(:,:,:,:,1);

m_air = 4.81e-26;
kB    = 1.38e-23;

Nv = 200;
v_rms_max = max( sqrt(kB * TemperatureDim(:) / m_air) );
vmax = 6 * v_rms_max;

v_grid = linspace(-vmax, vmax, Nv);
dv = v_grid(2) - v_grid(1);



tau = 0.01;
dt  = 0.01;
dx  = 1.0;

f = computeBGKf(TemperatureDim, Ux, Rho, tau, v_grid, dt, dx);
size(f)

v3 = reshape(v_grid, [1, Nv, 1]);

rho_from_f = squeeze( sum(f, 2) ) * dv;

rhou_from_f = squeeze( sum(v3 .* f, 2) ) * dv;
u_from_f = rhou_from_f ./ rho_from_f;

q_from_f = 0.5 * m_air * squeeze( sum((v3.^3) .* f, 2) ) * dv;

disp('Moment sizes:')
size(rho_from_f)
size(q_from_f)

rho_input = squeeze(Rho(:,1,1,:));
rho_err = rho_from_f - rho_input;

fprintf('max |rho - rho_input| = %.3e\n', max(abs(rho_err(:))));

x_index = 50;
time = (0:Nt-1) * dt;

figure;
plot(time, rho_input(x_index,:), 'k', ...
    time, rho_from_f(x_index,:), '--r');
legend('Input \rho', 'Recovered \rho');
xlabel('time'); ylabel('\rho');

figure;
plot(time, rho_err(x_index,:));
xlabel('time'); ylabel('\rho error');

ix = 50;
t  = 1000;

rho_loc = Rho(ix,1,1,t);
u_loc   = Ux(ix,1,1,t);
T_loc   = TemperatureDim(ix,1,1,t);



rho_mass = rho_from_f * m_air;

T0 = squeeze(TemperatureDim(:,1,1,1));

k_air  = 0.025;
cv_air = 718;



T_hist= solveHeat1D(T0, rho_mass, q_from_f, k_air, cv_air, dx, dt');
T_err = T_hist - T_ref;
rel_T_err = T_err ./ max(abs(T_ref), 1);
figure;
plot(time, rel_T_err(x_index,:));
xlabel('time');
ylabel('temperature error');
title('temperature error');

figure;
plot(time, T_ref(x_index,:), 'k', ...
    time, T_hist(x_index,:), '--r');
legend('Reference T', 'Recovered T )');
xlabel('time'); ylabel('Temperature');

ix = 50;
it = 1500;

f_slice = f(ix,:,it);

figure;
plot(v_grid, f_slice, 'LineWidth', 2);
xlabel('v');
ylabel('f(x,v,t)');
title(sprintf('Velocity distribution at x=%d, t=%d', ix, it));
grid on;


x_index = 50;
y_index = 1;
z_index = 1;
Temperature = reshape(Temperature, [100, 1, 1, 2000]);
particle_velocity = reshape(particle_velocity, [100, 1, 1, 2000]);


TemperatureDim = 3*Temperature + T   ;

m_air = 4.81e-26;  % kg
kB = 1.38e-23;     % J/K

Ux = U(:,:,:,:,1);
F_L96 = particle_velocity;
F_L96 = F_L96 - mean(F_L96(:));
F_L96 = F_L96 / std(F_L96(:));
v_rms = sqrt(kB * TemperatureDim / m_air);
particle_velocityDim = Ux + F_L96 .* v_rms  ;
T_series   = squeeze(TemperatureDim(x_index, y_index, z_index, :));
P_series   = squeeze(P(x_index, y_index, z_index, :));
Rho_series = squeeze(Rho(x_index, y_index, z_index, :));
Ux_series  = squeeze(Ux(x_index, y_index, z_index, :));
v_series   = squeeze(particle_velocityDim(x_index, y_index, z_index, :));

m = size(T,4);
dt = 0.01;
time = (0:m-1) * dt;

figure;

plot(time, T_series);
xlabel('Time [s]');
ylabel('Temperature [K]');
title('Temperature Time Series');

figure;
plot(time, P_series);
xlabel('Time [s]');
ylabel('Pressure [Pa]');
title('Pressure Time Series');

figure;
plot(time, Rho_series);
xlabel('Time [s]');
ylabel('Density [kg/m^3]');
title('Density Time Series');

figure;
plot(time, Ux_series);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('X-Velocity Time Series');

sgtitle('U velocity at x=50');


figure;
plot(time, v_series);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Partilce velocity Time Series');
