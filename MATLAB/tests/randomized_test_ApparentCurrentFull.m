function randomized_test_ApparentCurrentFull(N)
    for k = 1:N
        n = 2; m = 2; % keep small for simplicity
        percent_more_food = rand(n,m,4) - 0.5; % values in [-0.5,0.5]
        current = rand(n+1,m+1,2)*4 - 2;       % values in [-2,2]
        fish_speed = randi([1 5]);

        result = ApparentCurrentFull(percent_more_food, current, fish_speed);

        % Manually compute expected result
        expected = zeros(n,m,4);
        for i=1:n
            for j=1:m
                Ul = current(i,j,1);
                Vb = current(i,j,2);
                Ur = current(i+1,j,1);
                Vt = current(i,j+1,2);

                per_more_food = reshape(percent_more_food(i,j,:), 4, 1);

                unit_vector = zeros(size(per_more_food));
                unit_vector(per_more_food > 0) = 1;
                fourWayFishSpeed = fish_speed * unit_vector;

                dir_swim_speed_single = [ ...
                    max(fourWayFishSpeed(1) + Vt, 0), ...
                    max(fourWayFishSpeed(2) - Vb, 0), ...
                    max(fourWayFishSpeed(3) + Ur, 0), ...
                    max(fourWayFishSpeed(4) - Ul, 0) ];

                expected(i,j,:) = dir_swim_speed_single;
            end
        end

        if ~isequal(result, expected)
            error('Randomized test %d failed.\nExpected %s\nGot %s', ...
                k, mat2str(expected), mat2str(result));
        end
    end
end