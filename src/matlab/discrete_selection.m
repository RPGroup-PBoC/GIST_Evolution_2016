function p = discrete_selection(po, s, n_gen)
%discrete_selection loops through the number of generations n_gen and given
%an initial population proportion po and a selection coefficient it
%simulates consecutive generations to track how the proportion of both
%populations evolves over time.
%
% Parameter
% ---------
% po: float (0, 1].
%   Initial proportion of the "most fit population". Must be between 0 and
%   1
% s: float [0, 1].
%   Selection coefficient. Must be between 0 and 1
% n_gen: int.
%   Number of generations to run the simulation.
p = zeros(n_gen + 1, 1); % initialize array to save proportions
p(1) = po; % add the initial proportion
for i=1:n_gen
    p_i = p(i);
    q_i = 1 - p_i;
    p(i+1) = p_i + p_i * q_i * s / (1 - q_i * s);
end

end

