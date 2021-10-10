clear
clc
close all

%%
r = 0.015;
L = 1;
theta0 = 20;
tspan = [0,10];
iniCon = [0,deg2rad(theta0),0,0,1,0];
[t,t_step,xi,saturated_u,K_lqr,Psi,sum_state,Px,Py,M_pend] = Elastic_inv_pend_cart_ode(r,L,tspan,iniCon);

figure
subplot(4,1,1)
plot(t,xi(1,:))
ylabel('$x$ (m)', 'Interpreter','latex')
subplot(4,1,2)
plot(t,xi(4,:))
ylabel('$\dot{x}$ (m/s)', 'Interpreter','latex')
subplot(4,1,3)
plot(t,xi(2,:))
ylabel('$\theta$ (rad)', 'Interpreter','latex')
subplot(4,1,4)
plot(t,xi(5,:))
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter','latex')
xlabel('$t$ (sec)', 'Interpreter','latex')

figure
plot(t,saturated_u)
xlabel('$t$ (sec)', 'Interpreter','latex'); ylabel('$u$ (N)', 'Interpreter','latex')

figure
plot(t,Px(end,:))
xlabel('$t$ (sec)', 'Interpreter','latex'); ylabel('$P_x$ (m)', 'Interpreter','latex')

%% Animation
t_animation = t_step:t_step:tspan(end);
s = linspace(0,L,100);
bound = [-2 5 -0.125 L*2]; % boundary for display
xi_sim = [xi(1,:);xi(4,:);xi(2,:);xi(5,:)];
Animation_inverted_pendulum(xi_sim,Px,Py,s,t_animation,t_step,bound,'myVideo.avi');
