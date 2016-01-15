% -- Single experiment simulation -- %

% mutation rate
mut_rate = 3e-6;

% initial number of cells
No = 1000;

% number of generations
n_gen = 7;

initialize array of cells
n_mother = zeros(1, No);
for i=1:n_gen
    % generate random array
    n_daughter = rand(1, length(n_mother));
    % verify which cells are mutants
    n_daughter = n_daughter < mut_rate;
    % check if the mother cell was a mutant
    n_daughter = n_daughter | n_mother;
    n_mother = cat(2, n_mother, n_daughter);
end
n_mut = sum(n_mother);
%% 

% -- Many simulations -- %

% number of simulations
n_sim = 1000;

% initialize array to save outcomes of experiments
n_mut = zeros(1, n_sim);

for j=1:n_sim
    n_mother = zeros(1, No);
    for i=1:n_gen
        % generate random array
        n_daughter = rand(1, length(n_mother));
        % verify which cells are mutants
        n_daughter = n_daughter < mut_rate;
        % check if the mother cell was a mutant
        n_daughter = n_daughter | n_mother;
        %concatenate mother and daughter array
        n_mother = cat(2, n_mother, n_daughter);
    end
    n_mut(j) = sum(n_mother);
end

hist(n_mut, 100)

%%

% -- Include the possibility that BOTH daughter cells become mutants -- %

% initialize array to save outcomes of experiments
n_mut = zeros(1, n_sim);
for j=1:n_sim
    n_mother = zeros(1, No);
    for i=1:n_gen
        % generate random array
        n_daughter_1 = rand(1, length(n_mother));
        % verify which cells are mutants
        n_daughter_1 = n_daughter_1 < mut_rate;
        % check if the mother cell was a mutant
        n_daughter_1 = n_daughter_1 | n_mother;
        %daughter 2
        n_daughter_2 = rand(1, length(n_mother));
        % verify which cells are mutants
        n_daughter_2 = n_daughter_2 < mut_rate;
        % check if the mother cell was a mutant
        n_daughter_2 = n_daughter_2 | n_mother;
        %concatenate mother and daughter array
        n_mother = cat(2, n_daughter_1, n_daughter_2);
    end
    n_mut(j) = sum(n_mother);
end

hist(n_mut, 100)