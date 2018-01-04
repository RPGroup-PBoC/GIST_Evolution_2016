% simulate random mating between a population containing two alleles "1"
% and "2" of some gene

pop_size = 1000; % size of the population
n_gen = 200; % number of generations to run the simulation
n_loci = 2; % number of loci at which the gene appears

%probabilities of mating based on genotype
p_mate_homo1 = 0.4;
p_mate_hetero = 0.4;
p_mate_homo2 = 0.2;

% probabilities of different alleles in the first generation
p_1 = 0.5; % probability of allele "1" in the first generation
p_2 = 0.5; % probability of allele "2" in the first generation

pop_mat = zeros(n_gen,pop_size,n_loci);

% initialize the first generation of the population
for j=1:pop_size
    for l=1:n_loci
        r = rand();
        % randomly choose allele "1" or "2"
        if r < p_1
            pop_mat(1,j,l) = 1;
        else
            pop_mat(1,j,l) = 2;
        end
    end
end

% loop over subsequent generations of the simulation
for i=2:n_gen
    % loop over the individuals in the population
    for j=1:pop_size
        % randomly choose two individuals to reproduce
%         J_1 = randi(pop_size,1); % first randomly chosen to reproduce
%         J_2 = randi(pop_size,1); % second randomly chosen to reproduce
        p_reproduce = zeros(1,pop_size);
        % set p_reproduce to p_mate_homo1 everywhere we have homozygous 11
        p_reproduce(pop_mat(i-1,:,1)==1 & pop_mat(i-1,:,2)==1) = p_mate_homo1;
        % set p_reproduce to p_mate_hetero eveywhere we have heterozygous
        % 12 or 21
        p_reproduce(pop_mat(i-1,:,1)==1 & pop_mat(i-1,:,2)==2) = p_mate_hetero;
        p_reproduce(pop_mat(i-1,:,1)==2 & pop_mat(i-1,:,2)==1) = p_mate_hetero;
        % set p_reproduce to p_mate_homo2 everywhere we have homozygous 22
        p_reproduce(pop_mat(i-1,:,1)==2 & pop_mat(i-1,:,2)==2) = p_mate_homo2;
        
        % normalize p_reproduce
        p_reproduce = p_reproduce/sum(p_reproduce);
        
        J_1 = choose_int_from_distribution(p_reproduce); % first randomly chosen to reproduce
        J_2 = choose_int_from_distribution(p_reproduce); % second randomly chosen to reproduce
        
        % now choose one allele with equal probability from the first
        % parent
        r = rand();
        if r < 0.5
            pop_mat(i,j,1) = pop_mat(i-1,J_1,1);
        else
            pop_mat(i,j,1) = pop_mat(i-1,J_1,2);
        end
        
        % now choose one allele with equal probability from the second
        % parent
        r = rand();
        if r < 0.5
            pop_mat(i,j,2) = pop_mat(i-1,J_2,1);
        else
            pop_mat(i,j,2) = pop_mat(i-1,J_2,2);
        end
    end
end


n_homo_vec = zeros(1,n_gen);
n_homo1_vec = zeros(1,n_gen);
n_hetero_vec = zeros(1,n_gen);
n_homo2_vec = zeros(1,n_gen);

% loop over all individuals in all generations, and record whether they are
% homozygous or heterozygous
%first, loop over the generations
for i=1:n_gen
    % next, loop over individuals
    for j=1:pop_size
        % homozygous allele "1"
        if pop_mat(i,j,1)==1 && pop_mat(i,j,2)==1
            n_homo_vec(i) = n_homo_vec(i) + 1;
            n_homo1_vec(i) = n_homo1_vec(i) + 1;
        % heterzygous alleles
        elseif pop_mat(i,j,1) ~= pop_mat(i,j,2)
            n_hetero_vec(i) = n_hetero_vec(i) + 1;
        % finally, homozygous allele "2"
        else
            n_homo_vec(i) = n_homo_vec(i) + 1;
            n_homo2_vec(i) = n_homo2_vec(i) + 1;
        end
    end
end

% now plot the results of the simulation. First, we'll plot the percentage
% of the population that is homozygous and heterozygous as a function of
% time.
% we'll use subplot to lay the plots on top of one another
subplot(2,1,2)
plot(n_homo_vec/pop_size,'k')
hold on
plot(n_hetero_vec/pop_size,'g')
ylim([0 1])
hold off
xlabel('Generation number')
ylabel('Percent')
legend('percent homozygous', 'percent heterozygous')

% next, we'll use imagesc to make a visual map of how the population
% changes over time. We'll make use of a bit of hack to do this. Note that
% if we sum the values of the two different alleles, we'll get a unique
% number for each of the three cases: 1+1 = 2 (homozygous allele 1);
% 1+2 = 2+1 = 3 (heterozygous); 2+2 = 4. So if we use imagesc to plot the
% value of this sum for each individual, we'll get a unique color for each
% case. homozygous allele 1 will be blue, heterozygous will be green, and
% homozygous allele 2 will be red.
subplot(2,1,1)
imagesc(sort(sum(pop_mat,3),2)') 
xlabel('Generation number')
ylabel('Individual in population')
% in the above, transpose taken so that time runs left to right in the plot
% rather than top to bottom.