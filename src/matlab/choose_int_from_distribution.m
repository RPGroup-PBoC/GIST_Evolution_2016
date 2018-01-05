function i = choose_int_from_distribution(p)
% make sure p is normalized
p = p/sum(p);
% set up variables
c = cumsum(p);
chosen = false;
i = 1; % index of current choice
r = rand();
while ~chosen
    if r < c(i)
        chosen = true;
    else
        i = i+1;
    end
end
end