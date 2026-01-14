function [ dir_swim_speed_full ] = ApparentCurrentFull_happy( current, happiness, fish_speed, percent_more_food )
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
%  happinesss           : [n x m] binary array of hapiness
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
% DATE  : 30-12-2025
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
   
            % Is the fish happy where they are, used in swim logic later
            happy = happiness(i,j);

            % Is Hapiness a NaN? If so, quit now, no swimming
            if isnan( happy )
                dir_swim_speed_full(i,j,:) = [ 0 0 0 0 ];
                continue
            end

            Ul = current(i,j,1);
            Vb = current(i,j,2);
            Ur = current(i,j,1);
            Vt = current(i,j,2);
            
            dir_swim_speed_single = zeros(4,1);
            
            % Advection velocity, regardless of food
            
            % How fast to swim in each direction, based on percent more food

            per_more_food = reshape(percent_more_food(i,j,:), 4,1);

            % fourWayFishSpeed = mapSpeed( 'linear', per_more_food, fish_speed);
            
            % Need some logic that if the fish does not want to move out,
            % it stays where it is with some force?

            unit_vector = zeros(size(per_more_food));
            unit_vector( per_more_food > 0 ) = 1;
            fourWayFishSpeed = fish_speed * unit_vector;

            
            % Total speed, current PLUS swimming speed
            % This is a magnitude, which allows swimming against, BUT
            % also allows drift if current isoverpowering.
            % dir_swim_speed_single = [ abs( fourWayFishSpeed(1) + Vt ) ...
            %                           abs( fourWayFishSpeed(2) - Vb ) ...
            %                           abs( fourWayFishSpeed(3) + Ur ) ...
            %                           abs( fourWayFishSpeed(4) - Ul ) ];

            % dir_swim_speed = vel + fourWayFishSpeed;
            %% Swimming logic
            % Fish swim speed versus the current ( net gain or loss )
            vals = [  fish_speed - Vt ...
                      fish_speed + Vb ...
                      fish_speed - Ur ...
                      fish_speed + Ul ];
                        
            if happy
                % Swim speeds only where current is stronger than speed
                dir_swim_speed_single( vals < 0 ) = abs ( vals(vals<0) );
                
            else
                % dir_swim_speed_single = [ fourWayFishSpeed(1) + Vt ...
                %                           fourWayFishSpeed(2) + Vb ...
                %                           fourWayFishSpeed(3) + Ur ...
                %                           fourWayFishSpeed(4) + Ul ]
            
                % How fast are they swimming where they WANT to go
                dir_swim_speed_single = [ max( fourWayFishSpeed(1) + Vt, 0 ) ...
                                          max( fourWayFishSpeed(2) - Vb, 0 ) ...
                                          max( fourWayFishSpeed(3) + Ur, 0 ) ...
                                          max( fourWayFishSpeed(4) - Ul, 0 ) ];
                % Add the values where they are being pushed around by the current
                dir_swim_speed_single( vals <= 0 ) = abs ( vals(vals<=0) );
                
            end
            dir_swim_speed_full(i,j,:) = dir_swim_speed_single;

        end
    end
end