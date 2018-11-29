function fcollision = collisionCost(x, sig)
    global obstacles
    
    fcollision = 1;
    %% CHANGE THIS FOR 3D IMPLEMENTATION AND x SHOULD BE JUST THE STATES NOT VELOCITY
    for i = 1:size(obstacles, 2)
        xl = [obstacles{i}(1,1) obstacles{i}(2,1)];
        xu = [obstacles{i}(1,2) obstacles{i}(2,3)];
        sig = sig.^2;
        y = mvncdf(xl,xu,x',sig);
        fcollision = fcollision*y;
    end
end