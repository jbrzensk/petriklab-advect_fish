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

happy = 0;

if happy
    % Swim speeds onyl where current is stronger than speed
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


