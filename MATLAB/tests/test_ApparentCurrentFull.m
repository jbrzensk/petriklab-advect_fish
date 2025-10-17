function tests = test_ApparentCurrentFull
    % MATLAB unit test entry point (works with runtests)
    tests = functiontests(localfunctions);
end

%% ------------------------------------------------------------------------
function testFixedCase(testCase)
    % Fixed, hand-checkable example
    
    fish_speed = 2;
    
    % percent_more_food(i,j,:) = [Up, Down, Right, Left]
    percent_more_food = zeros(1,1,4);
    percent_more_food(1,1,1) = 0.5; % Up positive → swimming enabled
    percent_more_food(1,1,3) = 1.0; % Right positive → swimming enabled
    
    % current(i,j,:) = [u,v], pad to allow i+1,j+1 indexing
    current = zeros(2,2,2);
    current(1,1,:) = [1, 2];   % Ul = 1, Vb = 2
    current(2,1,:) = [3, 0];   % Ur = 3
    current(1,2,:) = [0, 4];   % Vt = 4
    
    result = ApparentCurrentFull(current, fish_speed, percent_more_food);
    
    % fourWayFishSpeed = [2, 0, 2, 0]
    % Up    = max(2 + Vt, 0) = max(2+4,0) = 6
    % Down  = max(0 - Vb, 0) = max(0-2,0) = 0
    % Right = max(2 + Ur, 0) = max(2+3,0) = 5
    % Left  = max(0 - Ul, 0) = max(0-1,0) = 0
    expected = reshape([6, 0, 5, 0], 1, 1, 4);
    
    verifyEqual(testCase, result, expected);
end

%% ------------------------------------------------------------------------
function testRandomizedCases(testCase)
    % Fuzz test with random inputs
    
    rng(42); % reproducible
    for k = 1:20
        n = 2; m = 2; % keep small for testing
        percent_more_food = randn(n,m,4);    % values ~N(0,1)
        current = randn(n+1,m+1,2);          % random currents
        fish_speed = randi([1 5]);
        
        result = ApparentCurrentFull(current, fish_speed, percent_more_food);
        
        % Compute expected manually (reference implementation)
        expected = zeros(n,m,4);
        for i=1:n
            for j=1:m
                Ul = current(i,j,1);
                Vb = current(i,j,2);
                Ur = current(i+1,j,1);
                Vt = current(i,j+1,2);

                per_more_food = reshape(percent_more_food(i,j,:),4,1);
                unit_vector = zeros(size(per_more_food));
                unit_vector(per_more_food > 0) = 1;
                fourWayFishSpeed = fish_speed * unit_vector;
            
                dir_swim_speed_single = [ abs( fourWayFishSpeed(1) + Vt ) ...
                                      abs( fourWayFishSpeed(2) - Vb ) ...
                                      abs( fourWayFishSpeed(3) + Ur ) ...
                                      abs( fourWayFishSpeed(4) - Ul ) ];
                % dir_swim_speed_single = [ ...
                %     max(fourWayFishSpeed(1) + Vt, 0), ...
                %     max(fourWayFishSpeed(2) - Vb, 0), ...
                %     max(fourWayFishSpeed(3) + Ur, 0), ...
                %     max(fourWayFishSpeed(4) - Ul, 0) ];
                
                expected(i,j,:) = dir_swim_speed_single;
            end
        end
        
        verifyEqual(testCase, result, expected, ...
            sprintf('Mismatch on randomized test %d', k));
    end
end
