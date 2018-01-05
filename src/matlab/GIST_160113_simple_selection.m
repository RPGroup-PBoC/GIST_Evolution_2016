% define number of generations
n_gen = 100;

% define selection coefficient
s_array = [0.1 0.05 0.005];

% define initial population fraction
po = 0.5;
qo = 1 - po;

% initialize array for population fraction
p_array = zeros(1, n_gen + 1);

% change value of first element of array
p_array(1) = po;

% run a for loop to update population fraction
for i=2:(n_gen + 1)
    p = p_array(i - 1);
    q = 1 - p;
    delta_p = p * q * s_array(1) / (1 - q * s_array(1));
    p_array(i) = p + delta_p;
end

%plot(1:(n_gen + 1), p_array)

% initialize population matrix
p_mat = zeros(length(s_array), n_gen + 1);

% update initial conditions
p_mat(:, 1) = po;

for i=1:length(s_array)
    for j=2:(n_gen + 1)
        p = p_mat(i, j - 1);
        q = 1 - p;
        delta_p = p * q * s_array(i) / (1 - q * s_array(i));
        p_mat(i, j) = p + delta_p;
    end
end

plot(1:n_gen+1, p_mat)
legend('0.1000', '0.05', '0.005')