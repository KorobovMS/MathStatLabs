function poisson = generate_poisson(lambda)
    p_of_i = @(k) lambda^k * exp(-lambda) / factorial(k);
    persistent probabilities;
    if (isempty(probabilities))
        probabilities = arrayfun(p_of_i, 0:1:100);
    end
    F = cumsum(probabilities);
    rnd = rand(1);
    poisson = find(F > rnd, 1) - 1;
end
