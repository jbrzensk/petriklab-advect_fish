function [ dir_swim_speed_full ] = ApparentCurrentFull_K( current, overcrowded, fish_speed, percent_more_food )
%% Calculates apparent current due to fish possibly swimming with 
% or against the current. Units should be the same of fish_speed
% and current
%
% INPUTS:
%  
%   current             : [m x n x 2] with fields 'u' and 'v' for velocity
%                         in x/y directions [m/s]
%                       : u = [ 1:m, 1:n, 1 ]
%                       : v = [ 1:m, 1:n, 2 ]
%  overcrowded          : [ n x m x 4] Overcrowded in each direction
%  fish_speed           : scalar, max swimming speed of fish [m/s]
%  percent_more_food    : [m x n x 4] array of food percentages in the four
%                         cardinal direction, up, down, right, left
%
% OUTPUT:
%   dir_swim_speed_full : [m x n x 4] array of swimming speeds, corrected
%                         for current, in the four cardinal direction 
%                         up, down, right, left.
%
% NOTES: 
%       - This array takes percent more food, in case we want to use a 
%         gradient for the swimming speed. We can adjust that here.
%
% AUTHOR: JARED BRZENSKI
% DATE  : 30-06-2025
% -------------------------------------------------------------------------
%% Core Apparent Current Full Code
    % Directions matrix
    % directions = [
    %             -1, 0;  % Up
    %              1, 0;  % Down
    %              0, 1;  % Right
    %              0, -1];% Left

    [n, m, ~] = size(percent_more_food);

    dir_swim_speed_full = zeros(n,m,4);

    for i=1:n
        for j=1:m
   
            % Check if in the ocean, Nan == not in ocean!
            if ( isnan( current(i,j,1) )) 
                dir_swim_speed_full(i,j,:) = [ 0 0 0 0 ];
                continue
            end

            Ul = current(i,j,1);
            Vb = current(i,j,2);
            Ur = current(i,j,1);
            Vt = current(i,j,2);
                        
            % val <= 0 means current pushing fish. Need to do that no matter what
            vals = [  fish_speed - Vt ...
                      fish_speed + Vb ...
                      fish_speed - Ur ...
                      fish_speed + Ul ];

            if (fish_speed <= 0 ) 
                dir_swim_speed_full(i,j,:) = vals;
                continue
            end

            % Dont care about directions I can fight against
            vals( vals>0 ) = 0;

            per_more_food = reshape(percent_more_food(i,j,:), 4,1);

            overcrowded_dir = reshape(overcrowded(i,j,:), 5, 1 );
            
            % Is the fish crowded where they are?
            crowded = overcrowded_dir(5);



            % Fish swim logic for overcrowded
            unit_vector = [ 0 0 0 0 ];
            
            if ( ~crowded )
                % Only try to swim towards more food
                unit_vector( per_more_food > 0 ) = 1;
            else
                % Swim in all directions, check which are crowded next
                unit_vector = [ 1 1 1 1];
            end
            
            unit_vector = unit_vector - overcrowded_dir(1:4)';

            % Only keep the directions we ACTUALLY want to go, not
            % negatives
            unit_vector = unit_vector == 1;
            
            % Swimming speed in direction I want to go, if no current
            want_to_swim = fish_speed * unit_vector;
            
            % How fast are they swimming where they WANT to go, if there is a
            %  current
            % If current faster than swim speed ( opposite ), returns 0.
            swim_with_the_current = unit_vector .* [  max( want_to_swim(1) + Vt, 0 ) ...
                                                      max( want_to_swim(2) - Vb, 0 ) ...
                                                      max( want_to_swim(3) + Ur, 0 ) ...
                                                      max( want_to_swim(4) - Ul, 0 ) ];

            % Add the values where they are being pushed around by the current
            % They cannot control this, so ,the zeros above might be getting pushed
            % have_to_swim = like_to_swim;
            have_to_swim = swim_with_the_current;
            have_to_swim(have_to_swim == 0) = abs(vals(have_to_swim == 0));
            

            dir_swim_speed_full(i,j,:) = have_to_swim;

        end
    end
end