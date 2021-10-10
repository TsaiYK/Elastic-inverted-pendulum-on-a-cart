function [t,t_step,xi,saturated_u,K_lqr,Psi,sum_state,Px,Py,M_pend] = Elastic_inv_pend_cart_ode...
    (r,L,tspan,iniCon)
% Independent parameters/variables
% r = 0.03;
density = 1.055e3;
para.E = 0.01e9;
para.u_max = 1000;
para.m = 0.1; % tip mass, unit: kg
para.M = 5; % cart mass, unit: kg
% para.L = 1;
para.L = L;
para.g = 9.81;
seg_x = 100;
theta0 = iniCon(2);
para.s = linspace(0,para.L,seg_x);
R_lqr = 0.0001;
% R = 0.001;
R = 0.0005;
% Q = diag([0.1,0.1,1,0.5,10,1]);
% Q = diag([0.1,0.1,1,1,10,1]);
Q = diag([0.1,0.01,5,0.5,10,1]);
Q_prime = 100*Q;

% Dependent parameters/variables
para.I = pi/4*r^4; % moment of area, unit: m^4
A = pi*r^2;
para.rho = density*A; % mass per length, unit: kg/m
M_pend = para.rho*para.L;
para.J = 1/3*M_pend*para.L^2;

[V,alpha] = elastic_analysis(theta0,seg_x,para);
para.alpha = alpha;
para.Vl = V(end);

%% feedback control
% K_lqr = lqr(A_dynamics,B_dynamics,Q,R);
A_approx = Linearization_A(0,0,theta0,0,0,0,para.alpha,V(end),diag(ones(1,6)*0.01),para);
B_approx = Linearization_B(0,0,theta0,0,0,0,0,para.alpha,V(end),0.01,para);
para.K = lqr(A_approx,B_approx,Q,R_lqr);
K_lqr = para.K;

[t, xi] = ode45(@(t,x) sys_cartpend(t,x,para), tspan, iniCon);
xi = xi';

% Redefine the state variables as [x,dx,theta,dtheta,T,dT]'.
state = [xi(1,:);xi(4,:);xi(2,:);xi(5,:);xi(3,:);xi(6,:)];
u = -K_lqr*state;
saturated_u = min(para.u_max, max(-para.u_max, u));

% v(s,t)=V(s)*T(t), where T(t)=xi(3,:)
V_x = V'*xi(3,:);

% cost func
t_step = mean(t(2:end)-t(1:end-1));
sum_state = sum(state(:,1:end-1).^2,2);
Psi = (sum(Q*sum_state)+sum(Q_prime*state(:,end).^2)+R*sum(u.^2))*t_step;

Px = -para.s'*sin(xi(2,:))-V_x.*cos(repmat(xi(2,:),100,1));
Py = para.s'*cos(xi(2,:))-V_x.*sin(repmat(xi(2,:),100,1));


    function dxi = sys_cartpend(t,xi,para)
    % Solve the dynamic equations by ode45
    % [Step 1]: define the generalized coordinates and control input
    % xi = [x,theta,T,dx,dtheta,dT]'
    x = xi(1);
    theta = xi(2);
    T = xi(3);
    dx = xi(4);
    dtheta = xi(5);
    dT = xi(6);
    
    % Full-state feedback controller u = -K*[x,dx,theta,dtheta,T,dT]
    u = -para.K*[x;dx;theta;dtheta;T;dT];
    saturated_u = min(para.u_max, max(-para.u_max, u));
    Vl = para.Vl;
    alpha = para.alpha;

    % [Step 2]: solve ddtheta and ddT by the linear system of equations
    a(1,1) = -(para.M+para.rho*para.L+para.m);
    a(2,1) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+para.m*Vl*T*sin(theta)+...
        para.rho*T*alpha(1)*sin(theta);
    a(3,1) = -para.m*cos(theta);

    b(1,1) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+...
        para.m*Vl*T*sin(theta)+para.rho*alpha(1)*T*sin(theta);
    b(2,1) = -para.J-para.m*para.L^2-para.rho*para.L^3/3-para.m*Vl^2*T^2-para.rho*alpha(2)*T^2;
    b(3,1) = -para.m*cos(theta);

    c(1,1) = -para.m*cos(theta);
    c(2,1) = 0;
    c(3,1) = -para.m*Vl;

    d(1,1) = para.m*Vl*dT*cos(theta)-2*para.m*dtheta*Vl*dT*sin(theta)-...
        para.m*dtheta^2*Vl*T*cos(theta)-(para.rho*para.L^2/2+para.m*para.L)*dtheta^2*sin(theta)-...
        2*para.rho*dtheta*alpha(1)*dT*sin(theta)-para.rho*alpha(1)*dtheta^2*T*cos(theta)-saturated_u;
    d(2,1) = 2*para.m*dtheta*Vl^2*T*dT-para.m*para.g*Vl*T*cos(theta)+...
        2*dtheta*para.rho*alpha(2)*T*dT-para.rho*para.g*T*cos(theta)-...
        (para.m*para.g*para.L+para.rho*para.g*para.L^2/2)*sin(theta);
    d(3,1) = -para.m*Vl*T*dtheta^2-para.m*para.g*sin(theta)-alpha(4)*T;

    % [Step 3]: Solve the system of algbraic equations (Cramer's law)
    Delta = det([a,b,c]);
    DeltaX = det([d,b,c]);
    DeltaY = det([a,d,c]);
    DeltaZ = det([a,b,d]);

    % [Step 4]: Obtain the derivative of xi(t): dxi(t)
    dxi(1,1) = dx;
    dxi(2,1) = dtheta;
    dxi(3,1) = dT;
    dxi(4,1) = DeltaX/Delta;
    dxi(5,1) = DeltaY/Delta;
    dxi(6,1) = DeltaZ/Delta;
