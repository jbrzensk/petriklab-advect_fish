function test_computeSwimSpeed()
    % Test parameters
    fish_speed = 2;

    % percent_more_food(i,j,:) corresponds to [Up, Down, Right, Left]
    % Let’s say only "Up" and "Right" have food (>0).
    percent_more_food = zeros(1,1,4);
    percent_more_food(1,1,1) = 0.5; % Up
    percent_more_food(1,1,3) = 1.0; % Right

    % current(i,j,:) = [U, V]
    % Pad with extra rows/cols since function calls i+1 and j+1
    current = zeros(2,2,2);
    current(1,1,:) = [1, 2];   % Ul=1, Vb=2
    current(2,1,:) = [3, 0];   % Ur=3
    current(1,2,:) = [0, 4];   % Vt=4

    % Run function
    result = ApparentCurrentFull(percent_more_food, current, fish_speed);

    % Expected:
    % fourWayFishSpeed = [2, 0, 2, 0]
    % Up    = max(2 + Vt, 0) = max(2+4,0) = 6
    % Down  = max(0 - Vb, 0) = max(0-2,0) = 0
    % Right = max(2 + Ur, 0) = max(2+3,0) = 5
    % Left  = max(0 - Ul, 0) = max(0-1,0) = 0
    expected = reshape([6, 0, 5, 0], 1, 1, 4);

    assert(isequal(result, expected), ...
        sprintf('Expected %s but got %s', mat2str(expected), mat2str(result)));

    disp('✅ test_computeSwimSpeed passed');
end
