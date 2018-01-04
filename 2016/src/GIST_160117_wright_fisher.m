% propagate the allele frequencies over generations using a wright fisher
% model

N = 16; % population size
n_gen = 20; % number of generations
P = zeros(2*N+1,2*N+1);

% initialize the transition matrix
for i=0:2*N
    for j=0:2*N
        P(j+1,i+1) = nchoosek(2*N,j)*((i/(2*N))^j)*((1-i/(2*N))^(2*N-j));
    end
end

X = zeros(2*N+1,n_gen);

% initialize the first vector
X(N+1,1) = 1.0;

% loop over the generations
for k=2:n_gen
    X(:,k) = P*X(:,k-1);
end

% imagesc(X)
% colorbar
bar3(X)