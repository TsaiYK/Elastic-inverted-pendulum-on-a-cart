function Animation_inverted_pendulum(xi_sim,Px,Py,X_vec,tspan,t_step,bound,filename)
% dim for the cart
Lx = 0.25;
Ly = 0.125;

y0 = 0;
h = figure;
set(gcf, 'Renderer', 'Painters');
% create the video writer with 1 fps
writerObj = VideoWriter(filename);
writerObj.FrameRate = 20;
% open the video writer
open(writerObj);

for i = 1:ceil(0.1/t_step):length(tspan)
    clf
    hold on
    grid on
    axis equal
%     set(gcf,'Position',[1,41,1920,963])
    
    axis(bound)
    title('Elastic Inverted Pendulum Animation')
    xlabel('x (m)','interpreter','latex')
    ylabel('y (m)','interpreter','latex')

    x0 = xi_sim(1,i);
%     xm = x0-X_vec*sin(xi_sim(2,i))-V_x(:,i)'*cos(xi_sim(2,i));
%     ym = y0+X_vec*cos(xi_sim(2,i))-V_x(:,i)'*sin(xi_sim(2,i));
    xm = x0+Px(:,i)';
    ym = y0+Py(:,i)';
    RecPlot(x0,y0,Lx,Ly)
    plot([x0,xm],[y0,ym],'k')
    for j = 2:length(X_vec)
        plot([xm(j-1),xm(j)],[ym(j-1),ym(j)],'k')
    end
    drawnow
    pause(0.01)
    % write the frames to the video
    F = getframe(gcf) ;
    writeVideo(writerObj, F);
end

% close the writer object
close(writerObj);
fprintf('Sucessfully generated the video\n')