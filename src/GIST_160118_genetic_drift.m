% define the population size
pop_size = 16; % as in Burin's exp

% define number of generations
n_gen = 30;

% initialize matrices to keep track of each fly's locus
locus_1 = zeros(n_gen, pop_size);
locus_2 = zeros(n_gen, pop_size);

% add alleles of generation 1.
locus_1(1, :) = ones(1, pop_size);
locus_2(1, :) = ones(1, pop_size) * 2;

% loop through generations selecting random flies to mate
for i=2:n_gen
    for j=1:pop_size
        % select random parents
        parents = randsample(1:pop_size, 2);
        % choose random alleles from each parent fly
        p_allele_1 = [locus_1(i-1, parents(1)) ...
            locus_2(i-1, parents(1))];
        p_allele_2 = [locus_1(i-1, parents(2)) ...
            locus_2(i-1, parents(2))];
        
        % select randomly the allele for the baby fly
        l_1 = randsample(p_allele_1, 1);
        l_2 = randsample(p_allele_2, 1);
        
        % save these alleles to keep track of them
        locus_1(i, j) = l_1;
        locus_2(i, j) = l_2;
    end
end

%%

pop_mat = sort((locus_1 + locus_2), 2);
subplot(3, 1, 1)
imagesc(pop_mat')

% plot genotypes frequencies
subplot(3, 1, 2)
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
subplot(3, 1, 3)
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

%%

% time to fixation as a function of population size

% define number of fly vials
n_vials = 107;

% define population size
pop_size_array = [4 8 16 32 64 128];

% initialize array to save # of gen to fix
n_gen_fix = zeros(length(pop_size_array), n_vials);

% loop through the number of population sizes
for p=1:length(pop_size_array)
    pop_size = pop_size_array(p)
    for i=1:n_vials
        % initialize arrays to keep track of each fly's loci
        locus_1 = ones(1, pop_size);
        locus_2 = ones(1, pop_size) * 2;
        
        % initialize boolean variable to use for while loop
        % this variable asks if one of the alleles was fixed or not
        fix_bool = all((locus_1 == 1 & locus_2 == 1) | ...
                       (locus_1 == 2 & locus_2 == 2));
        
        % initialize counting variable           
        n_gen = 1;
        
        % start while loop until fix_bool is true
        while fix_bool == false
            for j=1:pop_size
                % select random parents
                parents = randsample(1:pop_size, 2);
                % extract alleles from parents
                p_allele_1 = [locus_1(parents(1)) ...
                              locus_2(parents(1))];
                p_allele_2 = [locus_1(parents(2)) ...
                              locus_2(parents(2))];
                          
                % choose random allele from parents
                l_1 = randsample(p_allele_1, 1);
                l_2 = randsample(p_allele_2, 1);
                
                % save alleles to keep track of them
                locus_1(j) = l_1;
                locus_2(j) = l_2; 
            end
            fix_bool = all((locus_1 == 1 & locus_2 == 1) | ...
                       (locus_1 == 2 & locus_2 == 2));
            n_gen = n_gen + 1;
        end
        n_gen_fix(p, i) = n_gen;
    end
end
%%

plot(pop_size_array, mean(n_gen_fix, 2), '-o')
xlabel('Population size')
ylabel('Mean number of generations for fixation')