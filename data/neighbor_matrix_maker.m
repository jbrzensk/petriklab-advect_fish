%% Get all neighbors and save as mat mxnx4x2

close all; clear all; clc;

load('new_neighb_50x100.mat');
%load('index_neighborhood_100x50.mat')
%load('all_neighbors_2.mat');
%load('neighbors_360x200.mat');
%all_neighbors = neighborhood;
all_neighbors = new_neighborhood;

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

save('index_neighborhood_100x50_2.mat', 'neighborhood');