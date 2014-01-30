n = 10;
p = 0.1;
volume = 10000;
level_of_significance = 0.05;
gb = @(x) generate_binomial(n, p);

selection = arrayfun(gb, 1:1:volume);

%% Pearson criteria
histogram = hist(selection, n+1)';
empiric_possibility = histogram / volume;

%%
p_of_i = @(i) nchoosek(n, i) * p^i * (1-p)^(n-i);
theoretical_possibility = arrayfun(p_of_i, 0:1:n)';
    
%%
disp '--------------------------------'

summ = volume*sum( ((empiric_possibility-theoretical_possibility).^2) ./ theoretical_possibility )
chi = chi2inv(1 - level_of_significance, n)

if (summ <= chi)
    disp 'empiric distribution can be descibed as binomial';
else
    disp 'empiric distribution CANNOT be descibed as binomial';
end
