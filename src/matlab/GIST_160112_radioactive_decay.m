% Define the time
t = 0:0.01:3;

% Initialize vectors for element counts
N_K = zeros(1, length(t)); % Potassium
N_A = zeros(1, length(t)); % Argon

% Define the initial condition
N_Ko = 5E3;

% Define the decay constant
k = 5.81; % 10^-11 years^-1

% loop through time to calculate radioactive decay
for i=1:length(t)
   N_K(i) = N_Ko * exp(-k * t(i));
   N_A(i) = N_Ko * (1 - exp(-k * t(i)));
end

% Plot the element counts

plot(t, N_K, 'r')
hold on % keep the plot and add a new curve
plot(t, N_A, 'g')
plot(t, N_K + N_A, 'b')
% label the plot
xlabel('time 10^{-11} years^{-1}', 'FontSize', 20) % label the x axis
ylabel('Number of atoms', 'FontSize', 20) % label the y axis
legend('N_K', 'N_{Ar}', 'N_K + N_{Ar}')