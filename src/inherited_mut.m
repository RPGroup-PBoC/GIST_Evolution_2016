function num_mutants = inherited_mut(mut_rate, n_0, n_divisions)
% Function to simulate the inherited mutation hypothesis
% Parameter
% ---------
% mut_rate : float. 
%   probability of a basepair mutating
% n_0 : int. 
%   initial number of cells in the simulation
% n_divisions: int. 
%   number of simulated divisions
%
% Output
% ------
% num_mutants : int. 
%   Number of resistan mutants afteer n_divisions
cells = zeros(n_0, 1);
for i=1:n_divisions
    % duplicate the cells
    cells_temp = cells;
    % find cells that mutate generating random numbers between 0 and 1
    % those that are less than the mutation rate are considered resistant
    % mutants
    mutators = rand(length(cells_temp),1) < mut_rate; 
    % to account for cells whose mother cell was already a mutant we use
    % a simple boolean OR function. If the mother was mutant OR the
    % daughter turned to be a mutant the daughter must stay as mutant then
    cells_temp(mutators | cells_temp) = 1;
    % concatenate both arrays of daughter and mother cells
    cells = cat(1,cells,cells_temp);
end

num_mutants = sum(cells); % total number of mutants
end