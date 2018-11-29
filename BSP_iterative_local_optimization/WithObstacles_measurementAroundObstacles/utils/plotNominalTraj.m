% Plotting nominal trajectory
function plotNominalTraj(t0, x, u)
    global xf x0 T umax mapxmax mapymax
    
    set(gcf,'units','points','position',[10,10,700,700])
    
    subplot(323), hold on %PLOT control input ux
    plot(t0,u(1,:),'LineWidth',2);
    title(['Control (ux)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$u_x(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([-umax umax])
    hold on;

    subplot(324), hold on; %PLOT control input uy
    plot(t0,u(2,:),'LineWidth',2);
    title(['Control (uy)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$u_y(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([-umax umax])
    hold on;

    subplot(321), hold on %PLOT control input ux
    plot(t0,x(1,:),'LineWidth',2);
    title(['State (x)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$x(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([x0(1) xf(1)])
    hold on;

    subplot(322), hold on; %PLOT control input uy
    plot(t0,x(2,:),'LineWidth',2);
    title(['State (y)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$y(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([x0(2) xf(2)])
    hold on;

    subplot(3,2,[5,6]), hold on; %Phase Plot
    plot(x(1,:),x(2,:),'LineWidth',2);
    title(['Trajectory'], 'FontSize',24);
    grid on; xlabel('$x$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$y$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 mapxmax]); ylim([0 mapymax])
    hold on;
end