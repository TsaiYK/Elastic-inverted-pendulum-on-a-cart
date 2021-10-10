function A_approx = Linearization_A(x,dx,theta,dtheta,T,dT,alpha,Vl,delta_xi,para)
% The 'Linearization_A' function is to compute the approximate system matrix 
% whose form is linear.
for j = 1:6
    xi = [x;dx;theta;dtheta;T;dT];
%     u = -K_lqr*xi;
    u = 0;
    Aeq = zeros(3,3);
    Aeq(1,1) = -(para.M+para.rho*para.L+para.m);
    Aeq(1,2) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+...
        para.m*Vl*T*sin(theta)+para.rho*alpha(1)*T*sin(theta);
    Aeq(1,3) = -para.rho*alpha(1)*cos(theta);
    Aeq(2,1) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+para.m*Vl*T*sin(theta)+...
        para.rho*T*alpha(1)*sin(theta);
    Aeq(2,2) = -para.J-para.m*para.L^2-para.rho*para.L^3/3-para.m*Vl^2*T^2-para.rho*alpha(2)*T^2;
    Aeq(2,3) = -para.m*para.L*Vl-para.rho*alpha(3);
    Aeq(3,1) = -para.m*cos(theta);
    Aeq(3,3) = -para.m*Vl;
    Beq = zeros(3,1);
    Beq(1,1) = para.m*Vl*dT*cos(theta)-2*para.m*dtheta*Vl*dT*sin(theta)-...
        para.m*dtheta^2*Vl*T*cos(theta)-(para.rho*para.L^2/2+para.m*para.L)*dtheta^2*sin(theta)-...
        2*para.rho*dtheta*alpha(1)*dT*sin(theta)-para.rho*alpha(1)*dtheta^2*T*cos(theta)-u;
    Beq(2,1) = 2*para.m*dtheta*Vl^2*T*dT-para.m*para.g*Vl*T*cos(theta)+...
        2*dtheta*para.rho*alpha(2)*T*dT-para.rho*para.g*T*cos(theta)-...
        (para.m*para.g*para.L+para.rho*para.g*para.L^2/2)*sin(theta);
    Beq(3,1) = -para.m*Vl*T*dtheta^2-para.m*para.g*sin(theta)-alpha(4)*T;
    f(:,j) = Aeq\Beq;
    
    xi = xi+delta_xi(:,j);
    x = xi(1);
    dx = xi(2);
    theta = xi(3);
    dtheta = xi(4);
    T = xi(5);
    dT = xi(6);
    Aeq_delta = zeros(3,3);
    Aeq_delta(1,1) = -(para.M+para.rho*para.L+para.m);
    Aeq_delta(1,2) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+...
        para.m*Vl*T*sin(theta)+para.rho*alpha(1)*T*sin(theta);
    Aeq_delta(1,3) = -para.rho*alpha(1)*cos(theta);
    Aeq_delta(2,1) = -(para.m*para.L+para.rho*para.L^2/2)*cos(theta)+para.m*Vl*T*sin(theta)+...
        para.rho*T*alpha(1)*sin(theta);
    Aeq_delta(2,2) = -para.J-para.m*para.L^2-para.rho*para.L^3/3-para.m*Vl^2*T^2-para.rho*alpha(2)*T^2;
    Aeq_delta(2,3) = -para.m*para.L*Vl-para.rho*alpha(3);
    Aeq_delta(3,1) = -para.m*cos(theta);
    Aeq_delta(3,3) = -para.m*Vl;
    Beq_delta = zeros(3,1);
    Beq_delta(1,1) = para.m*Vl*dT*cos(theta)-2*para.m*dtheta*Vl*dT*sin(theta)-...
        para.m*dtheta^2*Vl*T*cos(theta)-(para.rho*para.L^2/2+para.m*para.L)*dtheta^2*sin(theta)-...
        2*para.rho*dtheta*alpha(1)*dT*sin(theta)-para.rho*alpha(1)*dtheta^2*T*cos(theta)-u;
    Beq_delta(2,1) = 2*para.m*dtheta*Vl^2*T*dT-para.m*para.g*Vl*T*cos(theta)+...
        2*dtheta*para.rho*alpha(2)*T*dT-para.rho*para.g*T*cos(theta)-...
        (para.m*para.g*para.L+para.rho*para.g*para.L^2/2)*sin(theta);
    Beq_delta(3,1) = -para.m*Vl*T*dtheta^2-para.m*para.g*sin(theta)-alpha(4)*T;
    f_delta(:,j) = Aeq_delta\Beq_delta;
    
    dfdxi(:,j) = (f_delta(:,j)-f(:,j))./delta_xi(j,j);
end

A_approx = [0,1,0,0,0,0;dfdxi(1,:);...
    0,0,0,1,0,0;dfdxi(2,:);...
    0,0,0,0,0,1;dfdxi(3,:)];
