% define population size
pop_size = 32;

% define number of generations
n_gen = 100;

% initialize transition matrix
P = zeros(2 * pop_size + 1, 2 * pop_size + 1);

% use for loops to build the matrix
for i=0:(2 * pop_size)
    for j=0:(2 * pop_size)
        P(j + 1, i + 1) = nchoosek(2 * pop_size, j) *...
            (i / (2 * pop_size))^j *...
            (1 - i / (2 * pop_size))^(2 * pop_size - j);
    end
end

% initialize vector for allele frequency
X = zeros(2 * pop_size + 1, n_gen);

% initilize the first vector
X(pop_size + 1, 1) = 1.0;

% apply the matrix on each generation
for k=2:n_gen
    X(:, k) = P * X(:, k - 1);
end

bar3(X)
ylabel('Number of alleles')
xlabel('Number of generations')
zlabel('Probability')