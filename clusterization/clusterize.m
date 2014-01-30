%% 
% Input: 
% -- dots - n-by-2 matrix of dots to clusterize
% -- means - k-by-2 matrix of centers of clusters
% -- dist - anonymous function handle @(x1, y1, x2, y2) of distance
%           between 2 dots
%%
function [clusters, means] = clusterize(dots, initial_means, dist)
dots_cardinality = size(dots, 1);
means_cardinality = size(initial_means, 1);

clusters = zeros(dots_cardinality, 1);
means = initial_means;
mc = 1:1:means_cardinality;
for i = 1:1:dots_cardinality    
    [~, ind] = min(dist(dots(i, 1), dots(i, 2), means(mc, 1), means(mc, 2)));
    clusters(i) = ind;
end

for i = 1:1:means_cardinality
    dots_in_cluster_i = dots(clusters == i, :);
    mean_x = mean(dots_in_cluster_i(:, 1));
    mean_y = mean(dots_in_cluster_i(:, 2));
    means(i, :) = [mean_x mean_y];
end
end
