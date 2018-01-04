% define population size
pop_size = 1000;

% define the number of generations to run the simulation
n_gen = 20;

% initialize matrices to keep track of each fly's alleles
locus_1 = zeros(n_gen, pop_size);
locus_2 = zeros(n_gen, pop_size);

% add alleles for 1st generation
locus_1(1, :) = ones(1, pop_size);
locus_2(1, :) = ones(1, pop_size) * 2;

% define weights
w11 = 1;
w12 = 2;
w22 = 0.5;

% STEP 1: loop through the generations
for i=2:n_gen
    % STEP 2: extract the genotypes of the previous generation
    p_gen = [locus_1(i - 1, :); locus_2(i - 1, :)];
    % STEP 3: compute the current weights
    p = sum(sum(p_gen == 1)) / (2 * pop_size);
    q = 1 - p;
    A11 = p^2 * w11;
    A12 = 2 * p * q * w12;
    A22 = q^2 * w22;
    w_bar = A11 + A12 + A22;
    % STEP 4: assign a weight to each fly from the previous generation
    p_fitness = zeros(1, pop_size);
    p_gen_sum = sum(p_gen, 1);
    p_fitness(p_gen_sum == 2) = A11 / w_bar;
    p_fitness(p_gen_sum == 3) = A12 / w_bar;
    p_fitness(p_gen_sum == 4) = A22 / w_bar;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
p = (sum(locus_1 == 1, 2) + sum(locus_2 == 1, 2)) /...
    (2 * pop_size);
delta_p = diff(p);
plot(delta_p, 'b')
xlabel('Number of generation')
ylabel('$\Delta p$', 'Interpreter', 'Latex')
xlim([1 n_gen])