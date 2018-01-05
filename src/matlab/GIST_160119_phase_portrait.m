% CASE 1
% define the weights
w11 = 10;
w12 = 4;
w22 = 2;

% define p vector
p = 0:0.01:1;

% initialize a vector to save ?p
delta_p = zeros(0, length(p));

% compute ?p
for i=1:length(p)
    P = p(i);
    Q = 1 - P;
    w_bar = P^2 * w11 + 2 * P * Q * w12 + Q^2 * w22;
    delta_p(i) = P / w_bar * (P * (w11 - w_bar) +...
                              Q * (w12 - w_bar));
end

plot(p, delta_p)
xlabel('$p$', 'Interpreter', 'Latex')
ylabel('$ \Delta p $', 'Interpreter', 'Latex')