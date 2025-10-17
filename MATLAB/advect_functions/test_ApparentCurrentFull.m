%% Test Apparent Current Full
% Tests apparent current due to fish possibly swimming with 
% or against the current. Units should be the same of fish_speed
% and current
%
% INPUTS:
%  
%   current             : [m x n x 2] with fields 'u' and 'v' for velocity
%                         in x/y directions [m/s]
%                       : u = [ 1:m, 1:n, 1 ]
%                       : v = [ 1:m, 1:n, 2 ]
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
% DATE  : 30-09-2025
% -------------------------------------------------------------------------
fish_speed = 1;
current = zeros(4,4,2);
percent_more_food = zeros(3,3,4);

% Check for no motion if no current and no more food
results = ApparentCurrentFull( current, fish_speed, percent_more_food );

assert( sum(results(:)) == 0, 'Zero current and food assertion fails' );

% Check for motion due to positive U
current(2,2,1) = 1;
results = ApparentCurrentFull( current, fish_speed, percent_more_food );

assert( results(1,2,3) == 1, 'Rightward current fails');
% Check for motion due to negative U
current(2,2,1) = -1;
results = ApparentCurrentFull( current, fish_speed, percent_more_food );

assert( results(2,2,4) == 1, 'Leftward current fails');

