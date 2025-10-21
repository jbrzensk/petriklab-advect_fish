function [ overcrowded_fish_full ] = OvercrowdedFishFourWay_K( Fish, neighborhood, K )
%% This finds fish in four directions and detects overcrowding
%
% INPUTS:
%   Fish         : [m x n] array of fish concentrations [unitless]
%   neighborhood : [n x m x 4 x 2] indices of neighboring cells
%   K            : carrying capacity of a cell
%
% OUTPUT:
%   more_fish_full : [ m x n x 5] binary ocercrowded fish (0 no, 1 yes )
%                                         in the four cardinal directions, 
%                                         up, down, left, and right, and
%                                         center
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
    %               0, -1, % Left
    %               0, 0]; % Curernt Cell

    [n, m] = size(Food);

    overcrowded_fish_full = zeros(n,m,5);

    for j=1:m
        for i=1:n
            
            current_val = Fish(i,j);
            % Reset food vector each loop
            crowded = [ 0 0 0 0 0];
            % Local Neighborhood
            neighbors = reshape(neighborhood(i,j,:,:), [4, 2]);

            % Loop through each direction ( find more food )
            for k = 1:4

                ni = neighbors(k, 1);
                nj = neighbors(k, 2);
                
                % Check bounds and compare food values
                if Fish(ni, nj) > K
                    crowded(k) = 1;
                end

            end
            % Cehck Self
            if ( Fish(i,j) > K )
                crowded(5) = 1;
            end

            overcrowded_fish_full(i,j,:) = crowded;

        end
    end
end
