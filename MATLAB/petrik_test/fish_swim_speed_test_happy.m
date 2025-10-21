% Swim Direction speed test

dir_swim_speed_single = [ 0 0 0 0 ];

swim_speed = 2;
fourWayFishSpeed = [ swim_speed 0 0 0 ];

Vt = 1;
Vb = -2;
Ur = 3;
Ul = -4;

% Fish swim speed versus the current ( net gain or loss )
vals = [  swim_speed - Vt ...
          swim_speed + Vb ...
          swim_speed - Ur ...
          swim_speed + Ul ];


            % Is the fish crowded where they are?
            crowded = overcrowded(i,j,5);

            Ul = current(i,j,1);
            Vb = current(i,j,2);
            Ur = current(i,j,1);
            Vt = current(i,j,2);
            
            dir_swim_speed_single = zeros(4,1);
            
            % Advection velocity, regardless of food
            
            % How fast to swim in each direction, based on percent more food

            per_more_food = reshape(percent_more_food(i,j,:), 4,1);

            overcrowded_dir = reshape(overcrowded(i,j,:), 5, 1 );
            
            % fourWayFishSpeed = mapSpeed( 'linear', per_more_food, fish_speed);
            

            % Need some logic that if the fish does not want to move out,
            % it stays where it is with some force?

            unit_vector = zeros(size(per_more_food));
            unit_vector( per_more_food > 0 ) = 1;

            % Subtract desire from overcrowded vector
            unit_vector = unit_vector - overcrowded_dir(1:4);
            % Only keep the directions we ACTUALLY want to go, not
            % negatives
            unit_vector = unit_vector == 1;

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
            
            happy = 0;
            
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



