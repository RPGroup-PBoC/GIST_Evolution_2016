% define the lambdas to use
lambda_array = [0.5 3 7];

% initialize values of k
k = 0:20;

% compute p(k)
%p = 3.^k * exp(-3) ./ factorial(k);


for i=1:length(lambda_array)
    p = lambda_array(i).^k * exp(-lambda_array(i)) ./ factorial(k);
    subplot(1, length(lambda_array), i)
    bar(k, p)
end
