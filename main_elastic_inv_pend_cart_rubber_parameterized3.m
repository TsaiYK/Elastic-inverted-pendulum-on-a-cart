% The script is to do the DoE of the performance of the elastic inverted
% pendulum on a cart given the parameters, L (length of the pendulum) and
% theta0 (initial anlge of the pendulum).
% An extra matric is added: mass of the pendulum, which is a physical-based
% attribute and it is time-independent
% The co-design problem becomes a multi-objective parametric co-design
% problem: fp(xp,para) = m_pend(r,l) & fc(xp,xc,para) = Psi(r,K,l,theta0).
% D.V.s: xp = r; xc = K;
% This framework is kind of a nested co-design.
% The model is referred to "A note on the vibration of a clamped-free beam 
% with a mass at the free end" by Laura et al. in 1974
% 
clear
clc
close all

% % Rubber
density = 1.055e3;

% r_vec = [0.005,0.01,0.05,0.1]; % radius of the cross section, unit: m
% r_vec = linspace(0.02,0.04,10); % radius of the cross section, unit: m
r_vec = linspace(0.01,0.05,20); % radius of the cross section, unit: m
% r_vec = linspace(0.01,0.05,50); % radius of the cross section, unit: m
% L_vec = 1.25; % length, unit: m
L_vec = linspace(0.9,1.1,3); % length, unit: m
% L_vec = 0.91:0.01:1.1; % length, unit: m
% theta0_vec = deg2rad(20);
theta0_vec = deg2rad([5,10,15,20]);
% theta0_vec = deg2rad(linspace(5,20,10));
% theta0_vec = deg2rad(5:1:20);

tspan = [0,10];
PlotOrNot = 'DontPlot';
seg_x = 100;
for k = 1:length(L_vec)
    for j = 1:length(theta0_vec)
        theta0 = theta0_vec(j);
        iniCon = [0,theta0,0,0,0,0];
        for i = 1:length(r_vec)
            A = pi*r_vec(i)^2;
            para.rho = density*A; % mass per length, unit: kg/m
            M_pend{1,k}(j,i) = para.rho*L_vec(k);
            [t{k,j,i},t_step{1,k}(j,i),xi{k,j,i},saturated_u{k,j,i},K_lqr{k,j,i},...
                Psi{1,k}(j,i),sum_state{k,j,i},Px{k,j,i},Py{k,j,i}] = ...
                Elastic_inv_pend_cart_ode(r_vec(i),L_vec(k),tspan,iniCon);
            if strcmp(PlotOrNot,'Plot')
                figure(1);
                subplot(4,1,1)
                plot(t{k,j,i},xi{k,j,i}(1,:),'LineWidth',2); hold on
                ylabel('$x(t)$ (m)','interpreter','latex')
                grid on
                
                subplot(4,1,2)
                plot(t{k,j,i},xi{k,j,i}(2,:),'LineWidth',2); hold on
                ylabel('$\dot{x}(t)$ (m/s)','interpreter','latex')
                grid on
                
                subplot(4,1,3)
                plot(t{k,j,i},rad2deg(xi{k,j,i}(3,:)),'LineWidth',2); hold on
                ylabel('$\theta(t)$ (deg)','interpreter','latex')
                grid on
                
                subplot(4,1,4)
                plot(t{k,j,i},xi{k,j,i}(4,:),'LineWidth',2); hold on
                xlabel('$t$ (sec)','interpreter','latex')
                ylabel('$\dot{\theta}(t)$ (rad/s)','interpreter','latex')
                grid on
                legendInfo{i} = ['r=' num2str(r_vec(i)) ' m'];
                
                figure(2);
                plot(t{k,j,i},saturated_u{k,j,i},'LineWidth',2); hold on
                xlabel('$t$ (sec)','interpreter','latex')
                ylabel('$u(t)$ (N$\cdot$m)','interpreter','latex')
                set(gcf, 'Renderer', 'Painters');
                legendInfo{i} = ['r=' num2str(r_vec(i)) ' m'];
            end
        end
%         [V_min(k,j),I_min(k,j)] = min(Psi{1,k}(j,:));
%         [V_max(k,j),I_max(k,j)] = max(Psi{1,k}(j,:));
%         I(k,j) = 20;

    end
end

if strcmp(PlotOrNot,'Plot')
    figure(1);
    legend(legendInfo)
    figure(2);
    legend(legendInfo)
end

%% Plot w.r.t radius
% figure
% plot(r_vec,Psi,'LineWidth',2);
% xlabel('$r$ (m)','interpreter','latex');ylabel('$\Psi$','interpreter','latex')

for k = 1:length(L_vec)
    for j = 1:length(theta0_vec)
        [V_min(k,j),I_min(k,j)] = min(Psi{1,k}(j,:),[],2);
        [V_max(k,j),I_max(k,j)] = max(Psi{1,k}(j,:));
        r_opt(k,j) = r_vec(I_min(k,j));
    end
end

%%
figure
[X,Y] = meshgrid(L_vec,theta0_vec);
X = X';
Y = Y';
surf(X,rad2deg(Y),r_opt)
xlabel('$l$ (m)','interpreter','latex')
ylabel('$\phi_0$ (deg)','interpreter','latex')
zlabel('$r^*_p$ (m)','interpreter','latex')


%% Plot responses
L_index = 1;
theta0_index = 1;
% r_index = [4,19];
% for theta0_index = 1:length(theta0_vec)
    for r_index = [I_min(L_index,theta0_index),I_max(L_index,theta0_index)]
    % for r_index = [1,I_min,length(r_vec)]
        figure((theta0_index-1)*2+1);
        sgtitle(['l=',num2str(L_vec(L_index)),' m, \phi_0=',...
            num2str(rad2deg(theta0_vec(theta0_index))),'^o'])
        subplot(4,1,1)
        plot(t{L_index,theta0_index,r_index},xi{L_index,theta0_index,r_index}(1,:),'LineWidth',2); hold on
        ylabel('$x(t)$ (m)','interpreter','latex')
        grid on

        subplot(4,1,2)
        plot(t{L_index,theta0_index,r_index},xi{L_index,theta0_index,r_index}(4,:),'LineWidth',2); hold on
        ylabel('$\dot{x}(t)$ (m/s)','interpreter','latex')
        grid on

        subplot(4,1,3)
        plot(t{L_index,theta0_index,r_index},rad2deg(xi{L_index,theta0_index,r_index}(2,:)),'LineWidth',2); hold on
        ylabel('$\theta(t)$ (deg)','interpreter','latex')
        grid on

        subplot(4,1,4)
        plot(t{L_index,theta0_index,r_index},xi{L_index,theta0_index,r_index}(5,:),'LineWidth',2); hold on
        xlabel('$t$ (sec)','interpreter','latex')
        ylabel('$\dot{\theta}(t)$ (rad/s)','interpreter','latex')
        grid on
        legend(['r_p=',num2str(r_vec(I_min(L_index,theta0_index)),'%.4f'),' m (min)'],...
            ['r_p=',num2str(r_vec(I_max(L_index,theta0_index)),'%.4f'),' m (max)'])
        set(gcf, 'Renderer', 'Painters');

        figure((theta0_index-1)*2+2);
        plot(t{L_index,theta0_index,r_index},saturated_u{L_index,theta0_index,r_index},'LineWidth',2); hold on
        xlabel('$t$ (sec)','interpreter','latex')
        ylabel('$u(t)$ (N$\cdot$m)','interpreter','latex')
        legend(['r_p=',num2str(r_vec(I_min(L_index,theta0_index)),'%.4f'),' m (min)'],...
            ['r_p=',num2str(r_vec(I_max(L_index,theta0_index)),'%.4f'),' m (max)'])
        title(['l=',num2str(L_vec(L_index)),' m, \phi_0=',...
            num2str(rad2deg(theta0_vec(theta0_index))),'^o'])
        set(gcf, 'Renderer', 'Painters');

    end
% end

% for r_index = [4,19]
%     figure(1); legend(['r_p=',num2str(r_vec(r_index)),' m'],['r_p=',num2str(r_vec(r_index)),' m'])
% end

%% Plot Psi w.r.t. rp
figure;
% k = 1;
% j = 12;
plot_num = 1;
for k = 1:length(L_vec)
    for j = 1:length(theta0_vec)
        subplot(length(L_vec),length(theta0_vec),plot_num)
%         subplot(4,4,plot_num)
        plot(r_vec,Psi{1,k}(j,:),'b.','Markersize',10); hold on
        plot(r_vec,Psi{1,k}(j,:),'b')
        xlabel('$r_p$ (m)','interpreter','latex')
        ylabel('$f$','interpreter','latex')
        plot_num = plot_num+1;
        title(['l=',num2str(L_vec(k)),' m, \phi_0=',...
            num2str(rad2deg(theta0_vec(j))),'^o'],...
            'Units', 'normalized'); % Set Title with correct Position
    end
end

%% Surface plots for different lengths
[X,Y] = meshgrid(r_vec,theta0_vec);
X = X';
Y = Y';
figure
for k = 1:length(L_vec)
    subplot(1,length(L_vec),k)
    surf(X,rad2deg(Y),Psi{1,k}')
    xlabel('$r_p$ (m)','interpreter','latex')
    ylabel('$\phi_0$ (deg)','interpreter','latex')
    zlabel('$f$','interpreter','latex')
    title(['l=',num2str(L_vec(k)),' m'])
    set(gcf,'Position',[40,432,1802,331])
end
    
    

%% Plot Psi/m_pend w.r.t. rp
% figure;
% plot_num = 1;
% for k = 1:length(L_vec)
%     for j = 1:length(theta0_vec)
%         subplot(length(L_vec),length(theta0_vec),plot_num)
%         plot(r_vec,f{1,k}(j,:))
%         xlabel('$r_p$ (m)','interpreter','latex')
%         ylabel('$\Psi/m_{pend}$','interpreter','latex')
%         plot_num = plot_num+1;
% %         axis([1 4 1.5 20])
%         title(['l=',num2str(L_vec(k)),' m, \phi_0=',...
%             num2str(rad2deg(theta0_vec(j))),'^o'],...
%             'Units', 'normalized'); % Set Title with correct Position
%         delete(legend)
%     end
% end

%% Animation
% i = I(k,j);
% L = L_vec(k);
% seg_x = 100;
% s = linspace(0,L,seg_x);
% elasticScale = 10;
% [t_test,t_step_test,xi_test,~,~,Psi,Px_test,Py_test,~] = Elastic_inv_pend_cart_ode_scalable...
%     (r_vec(i),L,tspan,iniCon,elasticScale);
% 
% 
% bound = [-2 5 -0.125 L*2]; % boundary for display
% xi_sim = [xi_test(1,:);xi_test(3,:);xi_test(2,:);xi_test(4,:)];
% Animation_inverted_pendulum3(xi_sim,Px_test,Py_test,s,t_test,t_step_test,bound,'myVideo.avi');

% figure;plot(x_V_rotate{1}(1,:),x_V_rotate{1}(2,:)); hold on
% plot(x_V_rotate{end}(1,:),x_V_rotate{end}(2,:))
% xlabel('$x$ (m)','interpreter','latex')
% ylabel('$y$ (m)','interpreter','latex')
% legend('initial','final')

%% Rotation mapping from the horizontal orientation to veritcal
% RotationMatrix = [cos(pi/2+theta),-sin(pi/2+theta);sin(pi/2+theta),cos(pi/2+theta)];
% x_V_rotate = RotationMatrix*[x;V];
% 
% figure
% % plot(x,X); hold on
% plot([0,0],[0,para.L],'k--'); hold on
% p1 = plot([0,para.L*cos(pi/2+theta)],[0,para.L*sin(pi/2+theta)],'k');
% p2 = plot(x_V_rotate(1,:),x_V_rotate(2,:),'b','LineWidth',2);
% axis equal
% % axis([0,L,X(1,end),0])
% set(gcf, 'Renderer', 'Painters');
% xlabel('$x$ (m)', 'Interpreter','latex')
% ylabel('$y$ (m)', 'Interpreter','latex')
% legend([p1,p2],'Rigid','Elastic')

