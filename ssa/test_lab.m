function test_lab()
% s = generate_time_series(1, 0, 0.1, 1, [1 2.1 2.9 4.8], 100);
% t = 1:numel(s);
t = 0:1:120;
% s = (t+0).*(1 + 0.05*(sin((1/2)*t).*cos((1/2)*t)+2*cos((1/2)*t)+3*sin(4*t)));
s = t + 5*rand(1, numel(t)).*sin(t);

figure;
hold on
grid on
plot(t, s, 'k')
plot(t, s, 'k*')
s = s(1:1:100);
for i = 100:1:120
    [~, forecast] = ssa(s(1:i), 50, 25, 12);
    plot(i, forecast, 'c*');
    s = cat(2, s, forecast);    
end

% [res, forecast] = ssa(s, 75, 25, 0);
% plot(t, res(1, :), 'r');
% plot(t, res(2, :), 'g');
% plot(t, res(3, :), 'b');

% [res2, ~] = ssa(res(1, :), 50, 25, 0);
% plot(t, res2(1, :), 'm');
% a = std(res(1, :))/std(t)
% 
% plot(numel(s)+1, forecast, 'c*');
% 
% for i = 1:1:50
%     s = cat(2, s, forecast);
%     [~, forecast] = ssa(s, 75, 25, 24);
%     %plot(numel(s)+1+i, forecast, 'c*');
% end
% figure
% plot(1:numel(s), s, 'k');

end

%% Generators
% y(k) = (a*t + b) + c*y(k-4) + d
%        --trend--   -season- -dev-
function series = generate_time_series(a, b, season_amp, dev_amp, init, volume)
series = zeros(1, volume);
init_size = numel(init);
indicies_to_copy = 1:init_size;
series(indicies_to_copy) = init(indicies_to_copy);

gen_rand_with_amp = @(amp) amp*(2*rand(1, 1)-1);

for i = init_size+1:volume
    series(i) = (a*i + b) + ...
        gen_rand_with_amp(season_amp)*series(i-4) + ...
        gen_rand_with_amp(dev_amp);
end
end

%% Singular spectrum analysis algorythm and its stages
% Description of variables is here:
% http://en.wikipedia.org/wiki/Singular_spectrum_analysis
function [decomposed_series, forecast] = ssa(x, L, m, r)
X = ssa_embedding(x, L);
[Xdecomp, forecast] = ssa_singular_value_decomposition(X, x, r);
Xgrouped = ssa_eigentriple_grouping(Xdecomp, m);
decomposed_series = ssa_averaging(Xgrouped);
end

function X = ssa_embedding(x, L)
N = numel(x);
K = N - L + 1;
X = zeros(L, K);
for i = 1:L
    X(i, :) = x(i:i+K-1);
end
end

function [Xdecomposition, forecast] = ssa_singular_value_decomposition(X, x, r)
S = X*X';
[vec, val] = eig(S);
to_sort = cat(2, diag(val, 0), vec');
sorted = flipud(sortrows(to_sort));
lambda = sorted(1:end, 1);
U = sorted(:, 2:end)';
d = sum(lambda > 0);
for i = 1:1:d
    V(:, i) = X'*U(:, i)/sqrt(lambda(i));
    Xdecomposition(i, :, :) = sqrt(lambda(i))*U(:, i)*V(:, i)';
end
if (r ~= 0)
    Vt = U(end, 1:r);
    V_star = U(1:end-1, 1:r);
    Q = x(numel(x) - numel(lambda) + 2 : numel(x))';
    forecast = Vt*V_star'*Q/(1 - Vt*Vt');
else
    forecast = 0;
end
end

function Xgrouped = ssa_eigentriple_grouping(Xdecomposition, m)
d = size(Xdecomposition, 1);
ind = 1:1:d;
n = hist(ind, m);
start_element = 1;
for i = 1:numel(n)
    ind_subset = ind(start_element:start_element+n(i)-1);
    start_element = start_element + n(i);
    if (numel(ind_subset) == 1)
        s = Xdecomposition(ind_subset, :, :);
    else
        s = sum(Xdecomposition(ind_subset, :, :));
    end
    Xgrouped(i, :, :) = reshape(s, size(s, 2), size(s, 3));
end
end

function decomposed_series = ssa_averaging(Xgrouped)
m = size(Xgrouped, 1);
for i = 1:1:m
    matrix = Xgrouped(i, :, :);
    height = size(matrix, 2);
    width = size(matrix, 3);
    matrix = reshape(matrix, height, width);
    matrix = fliplr(matrix);
    top_edge = width - 1;
    bottom_edge = height - 1;    
    for j = -bottom_edge:top_edge
        diagonal = diag(matrix, j);
        d_size = numel(diagonal);
        means = mean(mean(diagonal)) * ones(1, d_size);
        
        temp1 = diag(diagonal, j);
        temp1 = cat(1, temp1, zeros(height - size(temp1, 1), size(temp1, 2)));
        temp1 = cat(2, temp1, zeros(size(temp1, 1), width - size(temp1, 2)));
        temp1 = temp1(1:height, 1:width);
        
        temp2 = diag(means, j);
        temp2 = cat(1, temp2, zeros(height - size(temp2, 1), size(temp2, 2)));
        temp2 = cat(2, temp2, zeros(size(temp2, 1), width - size(temp2, 2)));
        temp2 = temp2(1:height, 1:width);
        
        matrix = matrix - temp1 + temp2;
    end
    matrix = fliplr(matrix);
    
    decomposed_series(i, :) = cat(2, matrix(:, 1)', matrix(end, 2:end));
end
end
