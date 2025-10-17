function neighbors_all = find_all_neighbors_by_index(ni, nj)
% Returns a struct of arrays, each ni x nj x 2, storing the [i,j] neighbors
% neighbors_all.north(i,j,:) = [row, col] of north neighbor
% neighbors_all.south(i,j,:) = ...
% neighbors_all.east(i,j,:)  = ...
% neighbors_all.west(i,j,:)  = ...
%
% NaN is used for neighbors outside the grid

directions = ["north", "south", "east", "west"];

% Initialize all fields with NaN arrays
for d = directions
    neighbors_all.(d) = NaN(ni, nj, 2);
end

% Fill in neighbors ( with periodic fillin)
for i = 1:ni
    for j = 1:nj
        % North
        if i < ni
            neighbors_all.north(i,j,:) = [i+1, j];
        else
            neighbors_all.north(i,j,:) = [1, j];
        end
        % South
        if i > 1
            neighbors_all.south(i,j,:) = [i-1, j];
        else
            neighbors_all.south(i,j,:) = [ni, j];
        end
        % East
        if j < nj
            neighbors_all.east(i,j,:) = [i, j+1];
        else
            neighbors_all.east(i,j,:) = [i, 1];
        end
        % West
        if j > 1
            neighbors_all.west(i,j,:) = [i, j-1];
        else
            neighbors_all.west(i,j,:) = [i, nj];
        end
    end
end

end
