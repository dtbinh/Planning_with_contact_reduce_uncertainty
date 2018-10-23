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
save('plan_outputs2.mat','armplan','armplanlength','armplanlength');
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
meanx = [];
meany = [];
% Create a movie structure to add frames to
mov(1:armplanlength) = struct('cdata', [], 'colormap', []);

% Create a video writer to use the writeVideo function
v = VideoWriter('ParticleRRTContacts2.avi');

% Make the video writer available for writing
open(v);
set(gcf,'units','points','position',[10,10,1000,1000])
for i = 1:armplanlength
    k = 1;
    for j = 1: size(armstart,2): numOfParticles*size(armstart,2)
        particlex(k) = armplan(i,j);
        particley(k) = armplan(i,j+1);
        k = k+1;
    end
    meanx = [meanx mean(particlex)];
    meany = [meany mean(particley)];
    
    scatter(particley,particlex,'filled');
    set(gca,'fontsize',18)
    text(armstart(2),armstart(1),'Start');
    text(armgoal(2),armgoal(1),'Goal');
    hold on;
    
    plot(meany, meanx, 'w-o', 'linewidth', 3);
    hold on;
    ax = gcf();
    mov(i) = getframe(ax);
%     [X,Map] = frame2im(mov(i));
% %     Write frame to the video writer "v"
    writeVideo(v,mov(i));
    pause(2);
end
%     movie(fig1,A,30,3,winsize)
    close(v);
end

%armplan