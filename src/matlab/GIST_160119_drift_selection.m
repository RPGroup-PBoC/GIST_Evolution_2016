% Define population size
pop_size = 1000; % as in Burin's experiment

% Define the number of generations to run the simulation
n_gen = 20; % also as in Burin's experiment

% initialize matrices to keep track of each fly's locus
locus_1 = zeros(n_gen, pop_size);
locus_2 = zeros(n_gen, pop_size);

% add alleles of generation 1.
locus_1(1, :) = ones(1, pop_size);
% for allele 2 we also use the function ones and simply multiply by 2.
locus_2(1, :) = ones(1, pop_size) * 2;

% define relative fitnesses
w11 = 1;
w12 = 2;
w22 = 2;

% STEP 1: loop through generations
for i=2:n_gen
    % STEP 2: extract genotypes of the previous generation
    p_gen = [locus_1(i-1, :); locus_2(i-1, :)];
    % STEP 3: compute the current weights
    p = sum(sum(p_gen == 1)) / (2 * pop_size);
    q = sum(sum(p_gen == 2)) / (2 * pop_size);
    A11 = p^2 * w11;
    A12 = 2 * p * q * w12;
    A22 = q^2 * w22;
    Z = A11 + A12 + A22;
    % STEP 4: assign a weight to each fly from the previous generation to
    % bias the parent selection.
    p_fitness = zeros(1, pop_size);
    p_fitness(sum(p_gen, 1) == 2) = A11 / Z;
    p_fitness(sum(p_gen, 1) == 3) = A12 / Z;
    p_fitness(sum(p_gen, 1) == 4) = A22 / Z;
    % STEP 5: loop through each of the newborn flies
    for j=1:pop_size
        % STEP 6: select random parents
        parents = randsample(1:pop_size, 2, true, p_fitness);
        % INTERMEDIATE STEP:
        % pick out the parents genotype from our locus_1 and locus_2 matrices
        p_allele_1 = [locus_1(i-1, parents(1)) ...
            locus_2(i-1, parents(1))];
        p_allele_2 = [locus_1(i-1, parents(2)) ...
            locus_2(i-1, parents(2))];
        % STEP 7: select randomly the allele for the baby fly
        l_1 = randsample(p_allele_1, 1);
        l_2 = randsample(p_allele_2, 1);
        % STEP 8: save these alleles to keep track of them
        locus_1(i, j) = l_1;
        locus_2(i, j) = l_2;
    end
end

%%
pop_mat = sort((locus_1 + locus_2), 2);
subplot(2, 2, 1)
imagesc(pop_mat')

% plot genotypes frequencies
subplot(2, 2, 2)
plot(sum(pop_mat == 2, 2) / pop_size, 'b')
hold on
plot(sum(pop_mat == 3, 2) / pop_size, 'g')
plot(sum(pop_mat == 4, 2) / pop_size, 'r')
legend('f_{11}', 'f_{12}', 'f_{22}')
xlabel('Number of generation')
ylabel('Frequency')
xlim([1 n_gen])
ylim([0 1])
hold off

% plot allele frequencies
subplot(2, 2, 3)
plot((sum(locus_1 == 1, 2) + sum(locus_2 == 1, 2)) /...
    (2 * pop_size), 'b')
hold on
plot((sum(locus_1 == 2, 2) + sum(locus_2 == 2, 2)) /...
    (2 * pop_size), 'r')
legend('f_1', 'f_2')
xlabel('Number of generation')
ylabel('Frequency')
xlim([1 n_gen])
ylim([0 1])
hold off

% plot ?p
subplot(2, 2, 4)
delta_p = diff((sum(locus_1 == 1, 2) + sum(locus_2 == 1, 2)) /...
    (2 * pop_size));
plot(delta_p, 'b')
xlabel('Number of generation')
ylabel('$\Delta p$', 'Interpreter', 'Latex')
xlim([1 n_gen])