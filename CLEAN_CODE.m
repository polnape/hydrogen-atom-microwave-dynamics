format long;
close all;
clear all;
clc
%%

%% HILLS REGIONS
% Parameters
% K = 0.1; % Parameter for the system
% h_values = [-1.7, -1.50155, -1.5]; % Energy levels to analyze

K = 0.0015749;
h_values = [-1.7, -1.421203, -1.35];
tol = 1e-13; % Convergence tolerance
max_iter = 1000; % Maximum iterations for Newton's method
delta_s = 1e-2; % Step size for pseudo-arcs
n = 2000; % Number of points per orbit
delta_y = 1e-4; % Increment for detecting sign changes in G
number = 50000;

% Compute equilibrium points
[x1, x2, h1, h2] = eq_and_energies(K);


x_range = linspace(-2, 2, 500);
y_range = linspace(-2, 2, 500);

tol = 1e-13; % Convergence tolerance
max_iter = 1000; % Maximum iterations for Newton's method
delta_s = 1e-2; % Step size for pseudo-arcs
n = 2000; % Number of points per orbit
delta_y = 1e-4; % Increment for detecting sign changes in G
number = 50000;



% Configuración del plot
figure;
for i = 1:length(h_values)
    h = h_values(i); % Energía actual

    % Inicializar la matriz de colores
    ColorMatrix = zeros(length(y_range), length(x_range));

    % Calcular el Hamiltoniano punto a punto
    for ix = 1:length(x_range)
        for iy = 1:length(y_range)
            energy = h + hamiltonian(x_range(ix), y_range(iy), K);
            if energy < 0
                ColorMatrix(iy, ix) = 1; % Color para valores permitidos
            else
                ColorMatrix(iy, ix) = -1; % Color para valores prohibidos
            end
        end
    end

    % Graficar las regiones
    subplot(1, 3, i);
    
    imagesc(x_range, y_range, ColorMatrix, [-1, 1]); % Limitar el rango de datos
    set(gca, 'YDir', 'normal'); % Invertir el eje Y para que sea correcto
    colormap([0 1 0; 1 0 0]); % Verde para permitido, rojo para prohibido
    title(sprintf('h = %.5f', h));
    xlabel('x');
    ylabel('y');
    axis equal;
    axis([-2 2 -2 2]); % Ajustar límites de los ejes
    grid on;
    hold on;
    plot(x1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(x2, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
end

h = -1.2;

% Inicializar la matriz de colores
ColorMatrix = zeros(length(y_range), length(x_range));

% Calcular el Hamiltoniano punto a punto
for ix = 1:length(x_range)
    for iy = 1:length(y_range)
        energy = h + hamiltonian(x_range(ix), y_range(iy), K);
        if energy < 0
            ColorMatrix(iy, ix) = 1; % Color para valores permitidos
        else
            ColorMatrix(iy, ix) = -1; % Color para valores prohibidos
        end
    end
end

% Crear la figura
figure;

% Graficar las regiones
imagesc(x_range, y_range, ColorMatrix, [-1, 1]); % Limitar el rango de datos
set(gca, 'YDir', 'normal'); % Invertir el eje Y para que sea correcto
colormap([0 1 0; 1 0 0]); % Verde para permitido, rojo para prohibido
title(sprintf('h = %.5f', h));
xlabel('x');
ylabel('y');
axis equal;
axis([-2 2 -2 2]); % Ajustar límites de los ejes
grid on;

% Añadir los puntos de equilibrio como crucetas negras
hold on;
plot(x1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2); % Puntos en (x1, 0)
plot(x2, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2); % Puntos en (x2, 0)
hold off;

%% EQ POINTS, EIGENVALUES AND EIGENVECTORS
% Define the parameter K
% K = 0.1; % Example value, you can modify it as needed
K =0.0015749;
% K = 0.028559865;

% Compute equilibrium points and energies
[x1, x2, h1, h2] = eq_and_energies(K);

% Display equilibrium points and their corresponding energies
disp('Equilibrium point x1:');
disp(x1);
disp('Energy at x1 (h1):');
disp(h1);

disp('Equilibrium point x2:');
disp(x2);
disp('Energy at x2 (h2):');
disp(h2);

% Calculate eigenvalues and eigenvectors for x1
[eigenvalues, eigenvectors] = equilibriumEigen(x1);
[eigenvalues2, eigenvectors2] = equilibriumEigen(x2);
% Display the results for x1
disp('Eigenvalues for x1:');
disp(eigenvalues);

disp('Eigenvectors for x1 (columns correspond to eigenvalues):');
disp(eigenvectors);

disp('Eigenvalues for x2:');
disp(eigenvalues2);

disp('Eigenvectors for x2 (columns correspond to eigenvalues):');
disp(eigenvectors2);

% Display the stable eigenvalue and eigenvector
disp('Stable eigenvalue (negative real part):');
disp(eigenvalues(1)); % First element corresponds to stable eigenvalue
disp('Stable eigenvector:');
disp(eigenvectors(:, 1)); % First column corresponds to stable eigenvector

% Display the unstable eigenvalue and eigenvector
disp('Unstable eigenvalue (positive real part):');
disp(eigenvalues(2)); % Second element corresponds to unstable eigenvalue
disp('Unstable eigenvector:');
disp(eigenvectors(:, 2)); % Second column corresponds to unstable eigenvector




v_unstable = eigenvectors(:, 2);
v_stable = eigenvectors(:, 1);
% Ensure eigenvectors have positive y-coordinate
if v_stable(2) > 0
    v_stable = -v_stable; % Flip stable eigenvector orientation
end

if v_unstable(2) > 0
    v_unstable = -v_unstable; % Flip unstable eigenvector orientation
end

%% MANIFOLDS
% Parámetros iniciales
s = 10e-4; % Escala para mover a lo largo del eigenvector inestable
i1 = [x1, 0, 0, x1] + s * v_unstable'; % Primer punto hacia el eigenvector inestable
i2 = [x1, 0, 0, x1] - s * v_unstable'; % Segundo punto en dirección opuesta
i3 = [x1, 0, 0, x1] +s*v_stable';
i4 = [x1, 0, 0, x1]-s*v_stable';
disp('Point along +eigenvector unstable:');
disp(i1);
disp('Point along -eigenvector unstable:');
disp(i2);

% Plotear los puntos
figure;
hold on;
scatter(i1(1), i1(2), 50, 'red', 'filled', 'DisplayName', 'Point +Eigenvector');
scatter(i2(1), i2(2), 50, 'yellow', 'filled', 'DisplayName', 'Point -Eigenvector');

% Etiquetas y formato
xlabel('x_1');
ylabel('x_2');
title('Initial Points Along the Unstable Eigenvector');
legend('show');
grid on;
hold off;


tol = 10e-2;
n_crossing = 2;
idir = 1;
[times1, sols1, orbit1] = poincare(@g, tol, K, i1, n_crossing, idir);
[times2, sols2, orbit2] = poincare(@g, tol, K, i2, n_crossing, idir);
[times3, sols3, orbit3] = poincare(@g, tol, K, i3, n_crossing, -idir);
[times4, sols4, orbit4] = poincare(@g, tol, K, i4, n_crossing, -idir);
% Plotear las órbitas

figure;
hold on;
plot(orbit1(:,1), orbit1(:,2), 'r-', 'DisplayName', 'Unstable branch (+)');
plot(orbit2(:,1), orbit2(:,2), 'magenta', 'DisplayName', 'Unstable branch (-)');
plot(orbit3(:,1), orbit3(:,2), 'b--', 'DisplayName', 'Stable branch (+)');
plot(orbit4(:,1), orbit4(:,2), 'g--', 'DisplayName', 'Stable branch (-)');

% Etiquetas y formato
xlabel('x_1');
ylabel('x_2');
title('Trajectories Along Stable and Unstable Eigenvectors');
legend('show');
grid on;
hold off;


disp(sols2);
%% LOOPS
% Parámetros iniciales
K_start = 0.028559865; % Valor inicial de K
K_step = 0.001;        % Incremento de K
K_end = 0.035;         % Valor final de K
tol = 10e-2;           % Tolerancia
n_crossing = 3;        % Número de cruces
idir = 1;              % Dirección del flujo
s = 10e-4;             % Escala para mover a lo largo del eigenvector inestable

% Inicializar figura
figure;
hold on;

% Iterar sobre valores de K
for K = K_start:K_step:K_end
    % Calcular puntos de equilibrio y energías
    [x1, x2, h1, h2] = eq_and_energies(K);

    % Calcular eigenvalores y eigenvectores para x1
    [eigenvalues, eigenvectors] = equilibriumEigen(x1);
    v_unstable = eigenvectors(:, 2); % Eigenvector inestable
    v_stable = eigenvectors(:, 1);   % Eigenvector estable

    % Asegurarse de que los eigenvectores tengan la orientación correcta
    if v_stable(2) > 0
        v_stable = -v_stable; % Cambiar orientación del estable
    end

    if v_unstable(2) > 0
        v_unstable = -v_unstable; % Cambiar orientación del inestable
    end

    % Punto inicial en la rama inestable (-)
    i2 = [x1, 0, 0, x1] - s * v_unstable';

    % Calcular órbita con Poincaré
    [~, ~, orbit2] = poincare(@g, tol, K, i2, n_crossing, idir);

    % Plotear la órbita
    plot(orbit2(:, 1), orbit2(:, 2), 'DisplayName', sprintf('K = %.3f', K));
end

% Etiquetas y formato
xlabel('x_1');
ylabel('x_2');
title('Trajectories for Varying K Values');
legend('show', 'Location', 'best');
grid on;
hold off;

%% PERIODIC ORBITS
xd = -0.5070094351999249; % Example x value
K = 0.0015749; % Given K
H = -1.7; % Given H
xri = 0.229690715908056;
xro = 2.155292149997760;

pyd = dy(xd, K, H);
pyri = dy(xri, K, H);
pyro = dy(xro, K, H);

disp(['p_y(1): ', num2str(pyd(1))]);
disp(['p_y(2): ', num2str(pyd(2))]);
disp(['p_y(1): ', num2str(pyri(1))]);
disp(['p_y(2): ', num2str(pyri(2))]);
disp(['p_y(1): ', num2str(pyro(1))]);
disp(['p_y(2): ', num2str(pyro(2))]);


%%
% Compute Poincare maps and orbits for each family
[PoincareMapTimes_d, sols_d, full_orbit_d] = poincare(@g, tol, K, [xd, 0, 0, pyd(2)], n_crossing, idir);
[PoincareMapTimes_ri, sols_ri, full_orbit_ri] = poincare(@g, tol, K, [xri, 0, 0, pyri(2)], n_crossing, idir);
[PoincareMapTimes_ro, sols_ro, full_orbit_ro] = poincare(@g, tol, K, [xro, 0, 0, pyro(2)], n_crossing, idir);
%%
% Plot full orbits for all families
figure;
hold on;
plot(full_orbit_d(:, 1), full_orbit_d(:, 2), 'g-', 'LineWidth', 2, 'DisplayName', 'Family d');
plot(full_orbit_ri(:, 1), full_orbit_ri(:, 2), 'b-', 'LineWidth', 2, 'DisplayName', 'Family ri');
plot(full_orbit_ro(:, 1), full_orbit_ro(:, 2), 'b--', 'LineWidth', 2, 'DisplayName', 'Family ro');
hold off;

% Add plot details
title('Full Orbits for Three Families (Up to 4 Crossings)', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
grid on;
legend show;
set(gca, 'FontSize', 12); % Make axis ticks larger
%% QUASIPERIODIC POINT
K = 0.0015749;
h = -1.7;
x = 0.1;
py = dy(x, K, H);
py_i = py(2);
tol = 10e-2;
[PoincareMapTimes_d, sols, full_orbit] = poincare(@g, tol, K, [x, 0, 0, py_i], 80, idir);
%%
figure;
hold on;
sols_filtered = sols(sols(:, 4) < 0, :);

plot(full_orbit(:, 1), full_orbit(:, 2), 'blue-', 'LineWidth', 0.5, 'DisplayName', 'Quasiperiodic orbit');
% scatter(sols_filtered(:, 1), sols_filtered(:, 2), 50, 'magenta', 'filled', 'DisplayName', 'Filtered crossings');
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
title('Quasi-periodic orbit, initial $(x,y) = (0.1,0)$', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

set(gca, 'FontSize', 14); % Make axis ticks larger
hold off;

%%
[t, sol] = ode45(@(t, state) f(t, state, K, h), tspan, init_cond);

%% LPO
h = -1.499;
K = 0.0015749;
% Compute equilibrium points and energies
[x1, x2, h1, h2] = eq_and_energies(K);

% Display equilibrium points and their corresponding energies
disp('Equilibrium point x1:');
disp(x1);
disp('Energy at x1 (h1):');
disp(h1);

disp('Equilibrium point x2:');
disp(x2);
disp('Energy at x2 (h2):');
disp(h2);

% Calculate eigenvalues and eigenvectors for x1
[eigenvalues, eigenvectors] = equilibriumEigen_imaginary(x1);
[eigenvalues2, eigenvectors2] = equilibriumEigen_imaginary(x2);

% Display the results for x1
disp('Eigenvalues for x1:');
disp(eigenvalues);

disp('Eigenvectors for x1 (columns correspond to eigenvalues):');
disp(eigenvectors);

disp('Eigenvalues for x2:');
disp(eigenvalues2);

disp('Eigenvectors for x2 (columns correspond to eigenvalues):');
disp(eigenvectors2);
%%
eps = 10^-3;
xinitial1 = x1+eps;
prime1 = dy(xinitial1, K, h);
ysign = +1;
tol = 10^-4;
max_counter = 10000;
der1 = prime1(2);
c = bisection_method(xinitial1, xinitial1, K, h, ysign, tol, max_counter);
%%
disp(c);
pyc = dy(c, K, h);
pyc = pyc(2);
[PoincareMapTimes_d, sols, full_orbit] = poincare(@g, tol, K, [x, 0, 0, pyc], 4, idir);
%% DISTANCES
K = 0.1;
h = -1.46536975;
[x1, x2, h1, h2] = eq_and_energies(K);
[eigenvalues, eigenvectors] = equilibriumEigen(x1);
[eigenvalues2, eigenvectors2] = equilibriumEigen(x2);
tol = 10e-2;
n_crossing = 2;
idir = 1;
s = 10e-4;
i2 = [x1, 0, 0, x1] - s * v_unstable';
[PoincareMapTimes_d, sols, full_orbit] = poincare_times(@g, tol, K, i2, 20, idir);
%%
[PoincareMapTimes, full_orbit] = poincare_times(@g, tol, K, i2, 300, idir);

% Plot r(t)
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column
r_t = sqrt(x_vals.^2 + y_vals.^2);

figure;
plot(times, r_t, 'LineWidth', 2);
title('Radial Distance r(t)', 'FontSize', 14);
xlabel('Time t', 'FontSize', 12);
ylabel('r(t)', 'FontSize', 12);
grid on;

%%
% Compute r(t) from full_orbit
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column

% Compute radial distance r(t)
r_t = sqrt(x_vals.^2 + y_vals.^2);
%%
% Plot r(t) with dashed lines
figure;
plot(times, r_t, '--', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Radial Distance $r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Time $t$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
legend('show', 'Interpreter', 'latex');% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels

% Plot the calculated orbit in the (x, y)-plane
figure;
plot(x_vals, y_vals, '-', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Trajectory in the $(x, y)$-Plane', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
legend('show', 'Interpreter', 'latex');
% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels
%%
K = 0.0015749;
per_radius = [xd, 0, 0, pyd(2)];
[PoincareMapTimes, full_orbit] = poincare_times(@g, tol, K, per_radius, 100, idir);


% Compute r(t) from full_orbit
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column

% Compute radial distance r(t)
r_t = sqrt(x_vals.^2 + y_vals.^2);
%%
% Plot r(t) with dashed lines
figure;
plot(times, r_t, '--', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Radial Distance $r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Time $t$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
ylim([0,1]);
legend('show', 'Interpreter', 'latex');% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels

% Plot the calculated orbit in the (x, y)-plane
figure;
plot(x_vals, y_vals, '-', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Trajectory in the $(x, y)$-Plane', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels

%%
% Initial condition
% Initial condition
init_cond = i2; % Replace with specific values
tspan = [0, 300]; % Extend the integration time range
K = 0.1; % Example parameter, adjust as needed
h = -1.7; % Energy level

% Options for ode45 with smaller time steps
tol = 1e-9; % Tighter tolerances for higher precision
options = odeset('RelTol', tol, 'AbsTol', tol, 'MaxStep', 1e-3); % Maximum step size for finer time resolution

% Integrate using ode45
[t, sol] = ode45(@(t, state) f(t, state, K, h), tspan, init_cond, options);

% Compute r(t) = sqrt(x(t)^2 + y(t)^2)
r_t = sqrt(sol(:, 1).^2 + sol(:, 2).^2);

% Plot r(t)
figure;
plot(t, r_t, 'LineWidth', 2);
title('Radial Distance r(t) with Fine Time Steps', 'FontSize', 14);
xlabel('Time t', 'FontSize', 12);
ylabel('r(t)', 'FontSize', 12);
grid on;
%% TWO BODY PROBLEM
K = 0;

h_values = [-1.7, -1.421203, -1.35];
tol = 1e-13; % Convergence tolerance
max_iter = 1000; % Maximum iterations for Newton's method
delta_s = 1e-2; % Step size for pseudo-arcs
n = 2000; % Number of points per orbit
delta_y = 1e-4; % Increment for detecting sign changes in G
number = 50000;

% Compute equilibrium points
[x1, x2, h1, h2] = eq_and_energies(K);
disp(h1);
disp(h2);

x_range = linspace(-2, 2, 500);
y_range = linspace(-2, 2, 500);

tol = 1e-13; % Convergence tolerance
max_iter = 1000; % Maximum iterations for Newton's method
delta_s = 1e-2; % Step size for pseudo-arcs
n = 2000; % Number of points per orbit
delta_y = 1e-4; % Increment for detecting sign changes in G
number = 50000;



% Configuración del plot
figure;
for i = 1:length(h_values)
    h = h_values(i); % Energía actual

    % Inicializar la matriz de colores
    ColorMatrix = zeros(length(y_range), length(x_range));

    % Calcular el Hamiltoniano punto a punto
    for ix = 1:length(x_range)
        for iy = 1:length(y_range)
            energy = h + hamiltonian(x_range(ix), y_range(iy), K);
            if energy < 0
                ColorMatrix(iy, ix) = 1; % Color para valores permitidos
            else
                ColorMatrix(iy, ix) = -1; % Color para valores prohibidos
            end
        end
    end

    % Graficar las regiones
    subplot(1, 3, i);
    
    imagesc(x_range, y_range, ColorMatrix, [-1, 1]); % Limitar el rango de datos
    set(gca, 'YDir', 'normal'); % Invertir el eje Y para que sea correcto
    colormap([0 1 0; 1 0 0]); % Verde para permitido, rojo para prohibido
    title(sprintf('h = %.5f', h));
    xlabel('x');
    ylabel('y');
    axis equal;
    axis([-2 2 -2 2]); % Ajustar límites de los ejes
    grid on;
    hold on;
    plot(x1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    plot(x2, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
end

%%
K = 0;
idir = +1;
x = 1.9;
py = dy(x, K, H);
py_i = py(2);
tol = 10e-2;
[PoincareMapTimes_d,full_orbit] = poincare(@g, tol, K, [x, 0, 0, py_i], 80, idir);


%%
figure;
hold on;
sols_filtered = sols(sols(:, 4) < 0, :);

plot(full_orbit(:, 1), full_orbit(:, 2), 'blue-', 'LineWidth', 0.5, 'DisplayName', 'Quasiperiodic orbit');
% scatter(sols_filtered(:, 1), sols_filtered(:, 2), 50, 'magenta', 'filled', 'DisplayName', 'Filtered crossings');
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
title('Quasi-periodic orbit, initial $(x,y) = (1.9,0)$', 'FontSize', 14, 'Interpreter', 'latex');
grid on;

set(gca, 'FontSize', 14); % Make axis ticks larger
hold off;

%% RADIAL DISTANCE DIF K
% K = 0.1;
K = 0.0015749;
% h = -1.46536975;
h = -1.4;
[x1, x2, h1, h2] = eq_and_energies(K);
disp(h1);
disp(h2);
[eigenvalues, eigenvectors] = equilibriumEigen(x1);
[eigenvalues2, eigenvectors2] = equilibriumEigen(x2);
tol = 10e-2;
n_crossing = 2;
idir = 1;
s = 10e-4;
% i2 = [x1, 0, 0, x1] - s * v_unstable';
x0 = 3.8;
py = dy(x0, K, h);
i2 = [x0, 0, 0, -py(2)];


[PoincareMapTimes, full_orbit] = poincare_times(@g, tol, K, i2, 1000, idir);
%%
% Plot r(t)
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column
r_t = sqrt(x_vals.^2 + y_vals.^2);

figure;
plot(times, r_t, 'LineWidth', 2);
title('Radial Distance r(t)', 'FontSize', 14);
xlabel('Time t', 'FontSize', 12);
ylabel('r(t)', 'FontSize', 12);
grid on;

%%
% Compute r(t) from full_orbit
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column

% Compute radial distance r(t)
r_t = sqrt(x_vals.^2 + y_vals.^2);
%%
% Plot r(t) with dashed lines
figure;
plot(times, r_t, '--', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Radial Distance $r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Time $t$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$r(t)$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
legend('show', 'Interpreter', 'latex');% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels

% Plot the calculated orbit in the (x, y)-plane
figure;
plot(x_vals, y_vals, '-', 'LineWidth', 1.5, 'DisplayName', '$W^{u}$ Branch');
title('Trajectory in the $(x, y)$-Plane', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
legend('show', 'Interpreter', 'latex');
% Adjust tick sizes
ax = gca; % Get current axes
ax.FontSize = 14; % Set font size for tick labels
%%
% Parámetros iniciales
% K = 0.0015749; % Valor dado de K
K = 0.1;
h = -1.7; % Valor de h para este plot
x0_range = linspace(1, 15, 500); % Valores de x0
theta_range = linspace(0, 2*pi, 500); % Valores de theta
[X0, Theta] = meshgrid(x0_range, theta_range); % Mallado para (x0, theta)

% Calcular velocidad inicial v según la relación
v_squared = 2 * (h + X0.^2 / 2 - 1 ./ X0 + K .* X0);
v_squared(v_squared < 0) = NaN; % Filtrar valores negativos para evitar complejos
V = sqrt(v_squared); % Velocidad inicial real

% Calcular energía osculante E_s^0
E_s = (V.^2 / 2) + (X0.^2 / 2) + X0 .* V .* sin(Theta) - (1 ./ X0);

% Figura para el plot
figure;
hold on;

% Región donde E_s^0 < 0 (magenta)
contourf(X0, Theta, real(E_s), [-Inf 0], 'm', 'DisplayName', 'E_s^0 < 0');

% Región donde E_s^0 > 0 (azul)
contourf(X0, Theta, real(E_s), [0 Inf], 'b', 'DisplayName', 'E_s^0 > 0');

% Formato del gráfico
title('Región Errática para h = -1.7');
xlabel('x_0');
ylabel('\theta');
legend show;
xlim([2,14]);
ylim([4, 5.5]);
grid on;
hold off;

%%
% Parámetros iniciales
h = -1.7; % Valor de h para este plot
x0_range = linspace(1, 15, 500); % Valores de x0
theta_range = linspace(0, 2*pi, 500); % Valores de theta
[X0, Theta] = meshgrid(x0_range, theta_range); % Mallado para (x0, theta)

% Valores de K para comparar
K_values = [0.001, 0.002, 0.005];
colors = [0.8 0 0.8; 0 0.8 0.8; 0.2 0.5 0.2]; % Colores lisos para las regiones de cada K

% Figura para el plot
figure;
hold on;

% Loop sobre los valores de K
for i = 1:length(K_values)
    K = K_values(i);

    % Calcular velocidad inicial v según la relación
    v_squared = 2 * (h + X0.^2 / 2 - 1 ./ X0 + K .* X0);
    v_squared(v_squared < 0) = NaN; % Filtrar valores negativos para evitar complejos
    V = sqrt(v_squared); % Velocidad inicial real

    % Calcular energía osculante E_s^0
    E_s = (V.^2 / 2) + (X0.^2 / 2) + X0 .* V .* sin(Theta) - (1 ./ X0);

    % Región donde E_s^0 < 0
    contourf(X0, Theta, real(E_s), [-Inf 0], 'LineStyle', 'none', 'FaceColor', colors(i, :), 'FaceAlpha', 0.7);
end

% Formato del gráfico
title('Región Errática para Diferentes Valores de K');
xlabel('x_0');
ylabel('\theta');
legend({'K = 0.001', 'K = 0.002', 'K = 0.005'}, 'Location', 'best');
xlim([2,14]);
ylim([4, 5.5]);
grid on;
hold off;

%%
K = 0.0015749;
h = -1.4;
[x1, x2, h1, h2] = eq_and_energies(K);
% [eigenvalues, eigenvectors] = equilibriumEigen(x1);
% [eigenvalues2, eigenvectors2] = equilibriumEigen(x2);
tol = 10e-2;
n_crossing = 2;
idir = 1;
s = 10e-4;
x0 = 8;
dyx0 = dy(x0, K, h);
i2 = [x0, 0, 0, -dyx0(2)];
disp(dyx0(2));

[PoincareMapTimes, full_orbit] = poincare_times(@g, tol, K, i2, 300, idir);

% Plot r(t)
times = full_orbit(:, 1); % Extract time column
x_vals = full_orbit(:, 2); % Extract x column
y_vals = full_orbit(:, 3); % Extract y column
r_t = sqrt(x_vals.^2 + y_vals.^2);

figure;
plot(times, r_t, 'LineWidth', 2);
title('Radial Distance r(t)', 'FontSize', 14);
xlabel('Time t', 'FontSize', 12);
ylabel('r(t)', 'FontSize', 12);
grid on;
%%
% Parámetros
K = 0.0015749;
h = -1.4;
x0 = 8;
theta = 4.85; % En radianes

% Cálculo de v
v = sqrt(2 * (h + 1/2 * x0^2 + 1/x0 - K * x0)); 

% Condiciones iniciales
xprima = cos(theta) * v;
yprima = sin(theta) * v;
initial = [x0, 0, xprima, -yprima];

% Opciones para ode45
tol = 1e-8; % Tolerancia
options = odeset('RelTol', tol, 'AbsTol', tol, 'MaxStep', 1e-3);

% Integración
[t, sol] = ode45(@(t, state) newf(t, state, K, h), [0, 100], initial, options);

% Graficar la trayectoria en el plano x-y
figure;
plot(sol(:,1), sol(:,2));
xlabel('x');
ylabel('y');
title('Trayectoria en el plano x-y');
grid on;


%%
function [x1, x2, h1, h2] = eq_and_energies(K)
    % Parameters
    tol = 1e-14; % Tolerance for Newton's method
    max_iter = 100; % Maximum number of iterations

    % Define f(x) and its derivative
    f1 = @(x) x^3 - K*x^2 + 1;
    df1 = @(x) 3*x^2 - 2*K*x;

    f2 = @(x) x^3 - K*x^2 - 1;
    df2 = @(x) 3*x^2 - 2*K*x;

    % Define Hamiltonian (energy function)
    H = @(x) 0.5 * x^2 + 1 / abs(x) - K * x;

    % Initial guesses based on the max expressions
    x1 = max(-1, -1/sqrt(K)); % Seed for x1 (L1)
    x2 = max(1, 2*K/3);       % Seed for x2 (L2)

    % Newton's method for x1
    for iter = 1:max_iter
        x1_new = x1 - f1(x1) / df1(x1);
        if abs(f1(x1_new)) < tol
            x1 = x1_new;
            break;
        end
        x1 = x1_new;
    end

    % Newton's method for x2
    for iter = 1:max_iter
        x2_new = x2 - f2(x2) / df2(x2);
        if abs(f2(x2_new)) < tol
            x2 = x2_new;
            break;
        end
        x2 = x2_new;
    end

    % Calculate energies
    h1 =-H(x1); % Energy at x1
    h2 = -H(x2); % Energy at x2

    % % Debugging Output
    % fprintf('Debugging Information:\n');
    % fprintf('x1 = %.14f, f1(x1) = %.14f\n', x1, f1(x1));
    % fprintf('x2 = %.14f, f2(x2) = %.14f\n', x2, f2(x2));
    % fprintf('h1 = %.14f (should be negative)\n', h1);
    % fprintf('h2 = %.14f (should be negative)\n', h2);
    % 
    % % Output results
    % fprintf('Equilibrium points and energies:\n');
    % fprintf('x1 = %.14f, h1 = %.14f\n', x1, h1);
    % fprintf('x2 = %.14f, h2 = %.14f\n', x2, h2);
end





function val = g(x)
    val = x(2);  % g(x) = velocity (second component of state vector)
end

function h = hamiltonian(x,y,K)
    r = (x^2+y^2)^(1/2);
    h = y^2/2+x^2/2+1/r-K*x;
 
end

function [eigenvalues, eigenvectors] = equilibriumEigen(xi)
    % equilibriumEigen - Computes and normalizes eigenvectors of the Jacobian
    % at the equilibrium point (xi, 0, 0, xi), and orders the eigenvalues and
    % eigenvectors such that:
    % - First is the stable eigenvalue and eigenvector (negative real part)
    % - Second is the unstable eigenvalue and eigenvector (positive real part)
    %
    % Input:
    %   xi - x-coordinate of the equilibrium point
    %
    % Output:
    %   eigenvalues - Vector with eigenvalues, first stable, second unstable
    %   eigenvectors - Matrix whose columns are the corresponding eigenvectors

    % Compute r = |xi|
    r = abs(xi);

    % Construct the Jacobian matrix
    J = [ 0,  1,         1,         0;
         -1,  0,         0,         1;
         (2 * xi^2) / r^5, 0,  0,         1;
         0, -1 / abs(xi)^3, -1,         0 ];

    % Compute eigenvalues and eigenvectors
    [raw_eigenvectors, D] = eig(J);
    raw_eigenvalues = diag(D); % Extract eigenvalues as a vector

    % Identify the indices for stable and unstable eigenvalues
    stable_index = find(real(raw_eigenvalues) < 0, 1); % First negative eigenvalue
    unstable_index = find(real(raw_eigenvalues) > 0, 1); % First positive eigenvalue

    % Extract and normalize the stable eigenvector
    stable_eigenvalue = raw_eigenvalues(stable_index);
    stable_eigenvector = raw_eigenvectors(:, stable_index) / norm(raw_eigenvectors(:, stable_index));

    % Extract and normalize the unstable eigenvector
    unstable_eigenvalue = raw_eigenvalues(unstable_index);
    unstable_eigenvector = raw_eigenvectors(:, unstable_index) / norm(raw_eigenvectors(:, unstable_index));

    % Construct the output lists
    eigenvalues = [stable_eigenvalue; unstable_eigenvalue];
    eigenvectors = [stable_eigenvector, unstable_eigenvector];
end


function [PoincareMapTimes, PoincareMapSols, full_orbit] = poincare(g, tol, K, x0, n_crossing, idir)
    % Maximum en uiterations allowed for refinement
    nmax = 100000;
    
    % Initialize outputs
    PoincareMapSols = zeros(n_crossing, 4);
    PoincareMapTimes = zeros(n_crossing, 1);
    full_orbit = []; % Full trajectory storage

    % ODE solver options
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

    % Time and step settings
    tau = 0;
    h = 0.01; % Initial step size
    tspan_fast = 0:1e-3:h;  % Larger steps for faster integration
    tspan_slow = 0:1e-4:0.001;  % Smaller steps for precise refinement

    % Plot initialization
    figure;
    hold on;
    grid on;
    xlabel('x');
    ylabel('y');
    title('Real-Time Orbit Plot in 2D');

    for ncross = 1:n_crossing
        if ncross ~= 1
            x0 = xk; % Update initial condition
        end

        % Initialize trajectory segment
        found = 0;
        orbit_segment = [];

        % Coarse integration until crossing is detected
        while ~found
            % Use larger steps before crossing
            [~, x_k1] = ode45(@(t, x) f(t, x, K, idir), [0, tspan_fast], x0, options);
            tau = tau + h;
            orbit_segment = [orbit_segment; x_k1]; % Append trajectory segment

            % Plot the segment
            plot(x_k1(:, 1), x_k1(:, 2), 'b');
            drawnow;

            % Check for crossing
            if g(x0) * g(x_k1(end, :)) < 0
                found = 1;

                % Exclude crossing point to avoid redundancy
                x_k1 = x_k1(1:end-1, :);
            end
            x0 = x_k1(end, :); % Update starting point for next step
        end

        % Refinement with smaller steps
        xk = x0; % Start from the detected point
        n = 0;
        while abs(g(xk)) > tol && n < 300
            % Compute time correction step
            delta_tau = -g(xk) / xk(4);
            tau = tau + idir * delta_tau;

            % Refine trajectory using smaller steps
            [~, x_k1] = ode45(@(t, x) f(t, x, K, idir), [0 abs(delta_tau)], xk, options);
            xk = x_k1(end, :);

            % Append refined points excluding duplicates
            orbit_segment = [orbit_segment; x_k1(1:end-1, :)];

            % Plot refined segment
            plot(x_k1(:, 1), x_k1(:, 2), 'b');
            drawnow;

            n = n + 1;
        end

        % Handle maximum iteration case
        if n >= nmax
            disp('Iteration surpassed');
            PoincareMapTimes(ncross, 1) = NaN;
            PoincareMapSols(ncross, :) = NaN;
        else
            % Store crossing information
            PoincareMapTimes(ncross, 1) = tau;
            PoincareMapSols(ncross, :) = xk;
        end

        % Append the trajectory segment to the full orbit
        full_orbit = [full_orbit; orbit_segment];
    end

    hold off;
end
function df = f(~, x, K, idir)
    df = zeros(4, 1);
    r = sqrt(x(1)^2 + x(2)^2);

    df(1) = x(3) + x(2);
    df(2) = x(4) - x(1);
    df(3) = x(4) - x(1) / r^3 - K;
    df(4) = -x(3) - x(2) / r^3;

    if idir == -1
        df = -df;
    end
end
% function py = dy(x, K, H)
%     py = -(2*H + x^2 + 2/abs(x) - 2*K*x)^(1/2);
% 
%     py = [py,py];
% end
function py = dy(x, K, H)
    % compute_py_branches - Computes the two possible values of p_y
    % based on the Hamiltonian equation.
    %
    % Inputs:
    %   x - Initial x value
    %   K - System parameter
    %   H - Energy level (Hamiltonian)
    %
    % Output:
    %   py - A vector containing the two possible p_y values:
    %        py(1): Positive root branch
    %        py(2): Negative root branch

    % Compute the discriminant
    discriminant = x^2 - 2 * (K * x - 1 / abs(x) - H);

    % Check if the discriminant is non-negative
    if discriminant < 0
        error('The discriminant is negative. No real solution exists for p_y.');
    end

    % Compute the two possible p_y values
    py_pos = x + sqrt(discriminant);  % Positive root branch
    py_neg = x - sqrt(discriminant);  % Negative root branch

    % Return the results as a vector
    py = [py_pos, py_neg];
end


% Bisection Method for Finding Periodic Orbits
function c = bisection_method(x1, x2, K, h, ysign, tol, max_counter)
    H = h;
    p_y1 = dy(x1, K, H);
    p_y1 = p_y1(2);
    p_y2 = dy(x2, K, H);
    p_y2 = p_y2(2);
    bool = 0;
    state1 = [x1, 0, 0, p_y1];
    state2 = [x2, 0, 0, p_y2];

    if F(x1, p_y1, K, H, tol) == 0
        c = x1;
        bool = 1;
    end
    if F(x2, p_y2, K, H, tol) == 0
        c = x2;
        bool = 1;
    end

    counter = 0;
    while abs(x2 - x1) > tol && counter < max_counter && bool == 0
        c = (x2 + x1) / 2;
        p_yc = dy(c, K, H);
        p_yc = p_yc(2);
        if F(c, p_yc, K, H, tol) == 0
            bool = 1;
        else
            if F(x1, p_y1, K, H, tol) * F(c, p_yc, K, H, tol) < 0
                x2 = c;
                p_y2 = p_yc;
            else
                x1 = c;
                p_y1 = p_yc;
            end
        end
        counter = counter + 1;
    end
    c = (x2 + x1) / 2;
end

function res = F(x, p_y, K, h, tol)
    x0 = [x, 0, 0, p_y];
    idir = +1;
    n_crossing = 1;
    [~, sols, ~] = poincare(@g, tol, K, x0, n_crossing, idir);
    res = sols(3); % Check third component (y') at Poincare crossing
end


% IMAGINARY PART ALSO EVEC, EVAL
function [eigenvalues, eigenvectors] = equilibriumEigen_imaginary(xi)
    % equilibriumEigen - Computes and returns all eigenvalues and eigenvectors
    % of the Jacobian at the equilibrium point (xi, 0, 0, xi).
    %
    % Input:
    %   xi - x-coordinate of the equilibrium point
    %
    % Output:
    %   eigenvalues - Vector with all eigenvalues
    %   eigenvectors - Matrix whose columns are the corresponding eigenvectors

    % Compute r = |xi|
    r = abs(xi);

    % Construct the Jacobian matrix
    J = [ 0,  1,         1,         0;
         -1,  0,         0,         1;
         (2 * xi^2) / r^5, 0,  0,         1;
         0, -1 / abs(xi)^3, -1,         0 ];

    % Compute eigenvalues and eigenvectors
    [eigenvectors, D] = eig(J);
    eigenvalues = diag(D); % Extract eigenvalues as a vector

    % Normalize eigenvectors
    for i = 1:size(eigenvectors, 2)
        eigenvectors(:, i) = eigenvectors(:, i) / norm(eigenvectors(:, i));
    end
end
    

function [PoincareMapTimes, full_orbit] = poincare_times(g, tol, K, x0, total_time, idir)
    % Computes the trajectory up to a specified total time without refinement.
    % Inputs:
    %   g          - Function defining the Poincaré section.
    %   tol        - Tolerance for checking crossings (still used for output).
    %   K          - System parameter.
    %   x0         - Initial condition (4D vector: [x, y, px, py]).
    %   total_time - Total simulation time.
    %   idir       - Direction of integration (+1 for forward, -1 for backward).
    %
    % Outputs:
    %   PoincareMapTimes - Times when the trajectory crosses the Poincaré section.
    %   full_orbit       - Complete trajectory data with time.

    % Initialize outputs
    PoincareMapTimes = [];
    full_orbit = [];

    % ODE solver options
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

    % Time and step settings
    tau = 0; % Accumulated time
    h = 10; % Initial step size
    tspan = [0, h]; % Time span for integration steps

    % Plot initialization
    figure;
    hold on;
    grid on;
    xlabel('x');
    ylabel('y');
    title('Real-Time Orbit Plot in 2D');

    % Simulate until total time is reached
    while tau < total_time
        % Integrate using ode45
        [t_segment, x_segment] = ode45(@(t, x)f(t, x, K, idir), tspan + tau, x0, options);

        % Update accumulated time and initial condition
        tau = t_segment(end);
        x0 = x_segment(end, :);

        % Store the trajectory segment
        full_orbit = [full_orbit; [t_segment, x_segment]];

        % Check for crossings with Poincaré section
        for i = 1:size(x_segment, 1) - 1
            if g(x_segment(i, :)) * g(x_segment(i + 1, :)) < 0
                % Linear interpolation for crossing time
                t_cross = t_segment(i) + (t_segment(i + 1) - t_segment(i)) * ...
                          (-g(x_segment(i, :))) / (g(x_segment(i + 1, :)) - g(x_segment(i, :)));
                PoincareMapTimes = [PoincareMapTimes; t_cross];
            end
        end

        % Plot the segment
        plot(x_segment(:, 1), x_segment(:, 2), 'b');
        drawnow;

        % Stop if total time is reached
        if tau >= total_time
            break;
        end
    end

    hold off;
end
function df = newf(~, x, K, idir)
    % Inicializamos df como un vector columna de ceros
    df = zeros(4, 1);
    
    % Cálculo de r
    r = sqrt(x(1)^2 + x(2)^2);
    
    % Asignamos las derivadas según el sistema reescrito
    df(1) = x(3); % x' = x(3)
    df(2) = x(4); % y' = x(4)
    df(3) = 2 * x(4) + x(1) - x(1) / r^3 - K; % x'' = 2y' + x - x/r^3 - K
    df(4) = -2 * x(3) - x(2) - x(2) / r^3;    % y'' = -2x' - y - y/r^3
    
    % Si idir es -1, invertimos el signo de df (para integración inversa)
    if idir == -1
        df = -df;
    end
end

