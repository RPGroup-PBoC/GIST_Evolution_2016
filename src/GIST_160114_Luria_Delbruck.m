% % probability that a particular cell will acquire the mutation
% mut_rate = 1e-5;
% % number of cells in the simulation
% n_cells = 128000;
% 
% % number of simulations
% n_simulations = 1000;
% % initialize array to save outcomes
% outcomes = zeros(n_simulations, 1);
% 
% for i=1:n_simulations
%     mut_cells = 0;
%    for j=1:n_cells
%        if rand(1, 1) < mut_rate
%            mut_cells = mut_cells + 1;
%        end
%    end
%    outcomes(i) = mut_cells;
% end

% mutation rate calibrated to give ? 1 mutant per plate.
mut_rate_inherit = 3e-6;
% initial number of cells
n_0 = 1000;
% number of cell divisions
n_divisions = 7;
% note that 128000 = 1000 * 2^7 if you were wondering why we chose that number

n_simulations = 1000;
outcomes_inherit = zeros(n_simulations, 1);

n_mother = zeros(1, n_0);
for j=1:n_simulations
    for i=1:n_divisions
        mutators = rand(1, length(n_mother)) < mut_rate_inherit;
        mutators = n_mother == 1 | mutators == 1;
        n_mother = cat(2, n_mother, mutators);
    end
    outcomes_inherit(j) = sum(n_mother);
end



