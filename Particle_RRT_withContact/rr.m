clc;
close all;
clear
% chdir /media/saumya/Data/Study/CMU/Courses/2ndSem/Planning_16-350/HW2/hw2code_16350_sp18/code
% mex planner_RRT.cpp
% mex planner_RRT.cpp
mex -v CXXFLAGS="$CXXFLAGS -std=c++11" -largeArrayDims particle_RRT.cpp tree.cpp
startQ = [30 30];
goalQ = [6 6];
runtest('map2.txt',startQ, goalQ);