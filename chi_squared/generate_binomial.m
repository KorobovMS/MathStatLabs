function binomial = generate_binomial(n, p)
    p_of_i = @(i) nchoosek(n, i) * p^i * (1-p)^(n-i);
    samples = 0:1:n;
    probabilities = arrayfun(p_of_i, samples);
    F = cumsum(probabilities);
    rnd = rand(1);
    binomial = find(F > rnd, 1) - 1;
end
