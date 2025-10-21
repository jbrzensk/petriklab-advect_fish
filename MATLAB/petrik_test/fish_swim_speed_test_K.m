% Swim Direction speed test
close all; clear all; clc;

dir_swim_speed_single = [ 0 0 0 0 ];

swim_speed = 2;
want_to_swim = [ swim_speed 0 0 0 ];
have_to_swim = [ 0 0 0 0 ];

per_more_food = [ 0 0 1 1 ];

% Is the fish crowded where they are?
overcrowded_dir = [ 0 0 0 0 1];
crowded = overcrowded_dir(5);


Vt = 1;
Vb = -2;
Ur = 3;
Ul = -4;

% Fish swim speed versus the current ( net gain or loss )
% val <= 0 means current pushing fish. Need to do that no matter what
vals = [  swim_speed - Vt ...
          swim_speed + Vb ...
          swim_speed - Ur ...
          swim_speed + Ul ];

% Dont care about directions I can fight against
vals( vals>0 ) = 0;

% Advection velocity, regardless of food

%unit_vector = zeros(size(per_more_food));

unit_vector = [ 0 0 0 0 ];

if ( ~crowded )
    % Only try to swim towards more food
    unit_vector( per_more_food > 0 ) = 1;
else
    % Swim in all directions, check which are crowded next
    unit_vector = [ 1 1 1 1];
end

unit_vector = unit_vector - overcrowded_dir(1:4);
% Only keep the directions we ACTUALLY want to go, not
% negatives
unit_vector = unit_vector == 1;

% Swimming speed in direction I want to go, if no current
want_to_swim = swim_speed * unit_vector;

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

%have_to_swim( vals < 0 ) = abs ( vals(vals<0) );
