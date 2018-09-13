function[numofmoves, caught] = runtest(mapfile, armstart, armgoal)

LINKLENGTH_CELLS=10;
envmap = load(mapfile);

close all;

%draw the environment
image(envmap, 'CDataMapping','scaled');
hold on;

%armplan should be a matrix of D by N 
%where D is the number of DOFs in the arm (length of armstart) and
%N is the number of steps in the plan 
[armplan, armplanlength, numOfParticles] = armplanner(envmap, armstart, armgoal); 

fprintf(1, 'plan of length %d was found\n', size(armplan,1));

%draw the plan
% midx = size(envmap,2)/2;
% x = zeros(length(armstart)+1,1);
% x(1) = midx;
% y = zeros(length(armstart)+1,1);
% fade = linspace(0.1,2,size(armplan,1));
% for i = 1:size(armplan)
%     for j = 1:length(armstart)
%         x(j+1) = x(j) + LINKLENGTH_CELLS*cos(armplan(i,j));
%         y(j+1) = y(j) + LINKLENGTH_CELLS*sin(armplan(i,j));
%     end;
%     plot(x,y, '-gs', 'linewidth', fade(i),...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);

%draw the plan
for i = 1:armplanlength
    k = 1;
    for j = 1: size(armstart,2): numOfParticles*size(armstart,2)
        particlex(k) = armplan(i,j);
        particley(k) = armplan(i,j+1);
        k = k+1;
    end
    scatter(particley,particlex,'filled');
    hold on;
    pause(2);
end
    
end

%armplan