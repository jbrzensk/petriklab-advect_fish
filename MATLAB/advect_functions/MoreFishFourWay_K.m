function [ more_fish_full ] = MoreFishFourWay_K( Fish, neighborhood, K )
%% This finds fish in four directions and detects overcrowding
%
% INPUTS:
%   Fish         : [m x n] array of fish concentrations [unitless]
%   neighborhood : [n x m x 4 x 2] indices of neighboring cells
%   K            : carrying capacity of a cell
%
% OUTPUT:
%   more_fish_full : [ m x n x 4] binary more fish (0 no, 1 yes )
%                                         in the four cardinal directions, 
%                                         up, down, left, and right
%
% Author: JARED BRZENSKI
% Date  : 20-10-2025
%
%  dir = fish MORE than current cell = [ U D R L ]
% -------------------------------------------------------------------------
%% Core percent more food code
    % Define neighbors relative to (i, j)
    % neighbors = [-1, 0;  % Up
    %               1, 0;  % Down
    %               0, 1;  % Right
    %               0, -1];% Left

    [n, m] = size(Food);

    more_fish_full = zeros(n,m,4);

    for j=1:m
        for i=1:n
            current_val = Fish(i,j);
            % Reset food vector each loop
            more_fish = [ 0 0 0 0 ];
            % Local Neighborhood
            % neighborhood = get_neighbors_from_struct(neighbors_all, i, j);
            % neighbors = [neighborhood.north; ...
            %              neighborhood.south; ...
            %              neighborhood.east; ...
            %              neighborhood.west];
            %neighbors = squeeze(neighborhood(i,j,:,:));
            neighbors = reshape(neighborhood(i,j,:,:), [4, 2]);

            % Loop through each direction ( find more food )
            for k = 1:4
                % ni = i + neighbors(k, 1); % Neighbor row
                % nj = j + neighbors(k, 2); % Neighbor column
                ni = neighbors(k, 1);
                nj = neighbors(k, 2);
                
                % Check bounds and compare food values
                %if ni >= 1 && ni <= n && nj >= 1 && nj <= m % Valid neighbor
                if Fish(ni, nj) > current_val
                    more_fish(k) = min( Food(ni,nj)/current_val, 2 ) - 1;
                end
                %end
            end
            more_fish_full(i,j,:) = more_fish;
        end
    end
end
