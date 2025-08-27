%% Get all neighbors and save as mat mxnx4x2

close all; clear all; clc;

load('all_neighbors_2.mat');
%load('neighbors_360x200.mat');

[m, n, ~] = size(all_neighbors.west);

neighborhood = zeros( m, n, 4, 2);

for ii=1:m
    for jj=1:n
        neighbors = get_neighbors_from_struct(all_neighbors, ii, jj);
        neighbor_vec = [neighbors.north; ...
                        neighbors.south; ...
                        neighbors.east; ...
                        neighbors.west];
        neighborhood( ii, jj, :, :) = neighbor_vec;
    end
end

save('all_neighbors_2_360x200.mat', 'neighborhood');