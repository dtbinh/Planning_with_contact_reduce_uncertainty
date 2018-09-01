function[armplan, armplanlength, numOfParticles] = armplanner(envmap, armstart, armgoal)


 %call the planner in C
 [armplan, armplanlength, numOfParticles] = particle_RRT(envmap, armstart, armgoal);

