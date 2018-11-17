function inside = checkIfInsideObstacle(x, y)
    global obstacles
    inside = 0;
    for i = 1:size(obstacles, 2)
        if x >= obstacles{i}(1,1) && x <= obstacles{i}(1,2) && y >= obstacles{i}(2,1) && y <= obstacles{i}(2,3)
            inside = 1;
        end
    end
        
end