function plotTrajectoryWithVariance(t0, b_new, u_new)

    global xf x0 T umax mapxmax mapymax nState N
    
    u = u_new;
    x = b_new(1:nState, :);
    sigmaVec = b_new(nState + 1:end, :);
%     set(gcf,'units','points','position',[10,10,700,700])
    
    subplot(423), hold on %PLOT control input ux
    plot(t0,u(1,:),'LineWidth',2);
    title(['Control (ux)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$u_x(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([-70 60])
    hold on;

    subplot(424), hold on; %PLOT control input uy
    plot(t0,u(2,:),'LineWidth',2);
    title(['Control (uy)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$u_y(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([-70 60])
    hold on;

    subplot(421), hold on %PLOT control input ux
    plot(t0,x(1,:),'LineWidth',2);
    title(['State (x)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$x(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([0 mapxmax/2])
    hold on;

    subplot(422), hold on; %PLOT control input uy
    plot(t0,x(2,:),'LineWidth',2);
    title(['State (y)'], 'FontSize',24);
    grid on; xlabel('$t$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$y(t)$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
    % axis square;
    xlim([0 T]); ylim([0 mapymax])
    hold on;

    subplot(4,2,[5,6,7,8]), hold on; %Phase Plot
    plot(x(1,:),x(2,:),'LineWidth',2);
    title(['Trajectory'], 'FontSize',24);
    grid on; xlabel('$x$', 'FontSize',24, 'Interpreter', 'latex'); ylabel('$y$', 'FontSize',24, 'Interpreter', 'latex');
    % xticks([-2*pi -pi 0 pi 2*pi]); xticklabels({'-2\pi','-\pi','0','\pi','2\pi'});
    set(gca,'fontsize',18)
%     axis square;
    xlim([0 mapxmax]); ylim([0 mapymax])
    hold on;
    for i = 1:9:N+1
        plot_gaussian_ellipsoid(x(:,i), vecTosigma(sigmaVec(:,i), nState))
        hold on;
    end
end