% define population size
pop_size = 16; % as in Burin's exp

% define the number of generations
n_gen = 50;

% initialize matrices to keep track of each fly's loci.
loci_1 = zeros(n_gen + 1, pop_size);
loci_2 = zeros(n_gen + 1, pop_size);

% define the genotype at the beginning of the experiment
loci_1(1, :) = ones(1, pop_size);
loci_2(1, :) = ones(1, pop_size) * 2;

% loop through generations selecting random flies to mate
for i=2:(n_gen + 1)
    for j=1:pop_size
        % select parents without replacement
        parents = datasample(1:16, 2, 'Replace', false);
        % extract the parents alleles
        p_allele_1 = [loci_1(i-1, parents(1)) loci_2(i-1, parents(1))];
        p_allele_2 = [loci_1(i-1, parents(2)) loci_2(i-1, parents(2))];

        % randomly choose alleles from the parents
        l_1 = datasample(p_allele_1, 1);
        l_2 = datasample(p_allele_2, 1);
        
        % save alleles to keep track of them
        loci_1(i, j) = l_1;
        loci_2(i, j) = l_2;
    end

end
sum(sum(pop_mat == 2))
%%

% we'll use imagesc to make a visual map of how the population changes over
% time. We'll make use of a bit of hack to do this. Note that if we sum the
% values of the two different alleles, we'll get a unique number for each
% of the three cases: 1+1 = 2 (homozygous allele 1); 1+2 = 2+1 = 3
% (heterozygous); 2+2 = 4. So if we use imagesc to plot the value of this
% sum for each individual, we'll get a unique color for each case.
pop_mat = sort((loci_1 + loci_2), 2);
subplot(3,1,1)
imagesc(pop_mat') 
xlabel('Generation number')
ylabel('Individual in population')
% in the above, transpose taken so that time runs left to right in the plot
% rather than top to bottom.

% finally, we will plot how the frequencies of the different alleles change
% over time
subplot(3,1,2)
% plot probability of homozygous 11
plot(sum(pop_mat==2, 2) / pop_size, 'b')
hold on
% plot probability of heterozygous 12 or 21
plot(sum(pop_mat == 3, 2) / pop_size, 'g')
% plot probability of homozygous 22
plot(sum(pop_mat == 4, 2) / pop_size, 'r')
xlabel('Generation number')
ylabel('Percent')
ylim([0 1])
xlim([1 n_gen])
legend('homozygous 11', 'heterozygous', 'homozygous 22','Location','Northwest')
hold off

% plot the allele probabilities
subplot(3,1,3)
plot((sum(loci_1==1, 2) + sum(loci_2==1, 2)) / (2 * pop_size), 'b');
hold on
plot((sum(loci_1==2, 2) + sum(loci_2==2, 2)) / (2 * pop_size), 'r');
xlabel('Generation number')
ylabel('Percent')
ylim([0 1])
xlim([1 n_gen])
legend('Allele 1','Allele 2')
hold off

%%

% Computing the fixation time
% to compute the fixation time we will use a while loop rather than a for
% loop.

% define the number of fly viles, i.e. experiments.
n_viles = 107;

% define population size
pop_size = 16; % as in Burin's exp

% initialize array to save fixation time
n_gens = zeros(1, n_viles);

for i=1:n_viles
    % initialize arrays to keep track of each fly's loci.
    loci_1 = ones(1, pop_size);
    loci_2 = ones(1, pop_size) * 2;
    fix_bool = all((loci_1 == 1 & loci_2 == 1) |...
                   (loci_1 == 2 & loci_2 == 2));
    n_gen = 1;
    while fix_bool == false;
        for j=1:pop_size
            % select parents without replacement
            parents = datasample(1:16, 2, 'Replace', false);
            % extract the parents alleles
            p_allele_1 = [loci_1(parents(1)) loci_2(parents(1))];
            p_allele_2 = [loci_1(parents(2)) loci_2(parents(2))];

            % randomly choose alleles from the parents
            l_1 = datasample(p_allele_1, 1);
            l_2 = datasample(p_allele_2, 1);

            % save alleles to keep track of them
            loci_1(j) = l_1;
            loci_2(j) = l_2;
        end
        fix_bool = all((loci_1 == 1 & loci_2 == 1) |...
                       (loci_1 == 2 & loci_2 == 2));
        n_gen = n_gen + 1;
    end
    n_gens(i) = n_gen;
    
end

hist(n_gens)
xlabel('Number of generations to fixation')
ylabel('Frequency')

%%

% define the number of fly viles, i.e. experiments.
n_viles = 107;

% define population size
pop_size_array = [4 8 16]; % as in Burin's exp

% initialize array to save fixation time
n_gens = zeros(length(pop_size_array), n_viles);

for p=1:length(pop_size_array)
    pop_size = pop_size_array(p)
        for i=1:n_viles
        % initialize arrays to keep track of each fly's loci.
        loci_1 = ones(1, pop_size);
        loci_2 = ones(1, pop_size) * 2;
        fix_bool = all((loci_1 == 1 & loci_2 == 1) |...
                       (loci_1 == 2 & loci_2 == 2));
        n_gen = 1;
        while fix_bool == false;
            for j=1:pop_size
                % select parents without replacement
                parents = datasample(1:pop_size, 2, 'Replace', false);
                % extract the parents alleles
                p_allele_1 = [loci_1(parents(1)) loci_2(parents(1))];
                p_allele_2 = [loci_1(parents(2)) loci_2(parents(2))];

                % randomly choose alleles from the parents
                l_1 = datasample(p_allele_1, 1);
                l_2 = datasample(p_allele_2, 1);

                % save alleles to keep track of them
                loci_1(j) = l_1;
                loci_2(j) = l_2;
            end
            fix_bool = all((loci_1 == 1 & loci_2 == 1) |...
                           (loci_1 == 2 & loci_2 == 2));
            n_gen = n_gen + 1;
        end
        n_gens(p, i) = n_gen;
    end
end

% plot the histogram for each

subplot(1, 2, 1)
hold on
for i=1:(size(n_gens, 1))    
    histogram(n_gens(i, :), 1:max_gen, 'facealpha', 0.3, ...
              'edgecolor', 'none');
end
legend('4', '8', '16')
xlabel('Number of generations to fixation')
hold off

subplot(1, 2, 2)
hold on
for i=1:(size(n_gens, 1))
    cdfplot(n_gens(i, :))
end
legend('4', '8', '16')
xlabel('Number of generations to fixation')
hold off