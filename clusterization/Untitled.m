x = [normrnd(2, 2.7, 100, 1); normrnd(20, 2.7, 100, 1); normrnd(2, 2.7, 100, 1); normrnd(10, 5.7, 100, 1)];
y = [normrnd(2, 2.7, 100, 1); normrnd(20, 2.7, 100, 1); normrnd(20, 2.7, 100, 1); normrnd(10, 5.7, 100, 1)];
dots = cat(2, x, y);
init_means = [5 5; 15 15; 5 15; 20 10];
dist = @(x1, y1, x2, y2) sqrt((x1 - x2).^2 + (y1 - y2).^2);
means = init_means;
for i = 1:1:10
    [clusters, means] = clusterize(dots, means, dist);



figure;
hold on;
for i = 1:1:400
    if (clusters(i) == 1)
        col = '.r';
    elseif (clusters(i) == 2)
        col = '.y';
    elseif (clusters(i) == 3)
        col = '.g';
    elseif (clusters(i) == 4)
        col = '.m';
    end
    plot(dots(i, 1), dots(i, 2), col);
end

plot(means(:, 1), means(:, 2), 'k*')
for i = 1:1:size(means, 1)
    tmp = dots(clusters == i, :);
    cv = convhull(tmp(:, 1), tmp(:, 2));
    plot(tmp(cv, 1), tmp(cv, 2), 'b')
end
end