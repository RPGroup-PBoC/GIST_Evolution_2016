% Define radius
r = 1:1000000; % meters

% Initialize vector to save energies
energy = zeros(1, length(r));

% Define constants
v = 10^3; % m / s
rho = 3000; % kg / m^3

% loop through r to calculate energy
for i=1:length(r)
    energy(i) = 2 * rho * r(i)^3 * v^2;
end

% Plot on a log-log scale
loglog(r, energy)
xlabel('radius [m]', 'FontSize', 20)
ylabel('energy [J]', 'FontSize', 20)