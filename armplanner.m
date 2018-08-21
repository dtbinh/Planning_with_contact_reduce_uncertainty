function[armplan] = armplanner(envmap, armstart, armgoal);


 %call the planner in C
 [armplan, armplanlength] = planner_RRT(envmap, armstart, armgoal);

