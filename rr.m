clc;
close all;
clear
% chdir /media/saumya/Data/Study/CMU/Courses/2ndSem/Planning_16-350/HW2/hw2code_16350_sp18/code
% mex planner_RRT.cpp
% mex planner_RRT.cpp
mex -v CXXFLAGS="$CXXFLAGS -std=c++11" -largeArrayDims planner_RRT.cpp tree.cpp
startQ = [pi/2 pi/4 pi/2 pi/4 pi/2];
goalQ = [pi/8 3*pi/4 pi 0.9*pi 1.5*pi];
runtest('map1.txt',startQ, goalQ);