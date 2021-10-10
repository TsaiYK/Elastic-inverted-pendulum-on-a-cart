function B_approx = Linearization_B(x,dx,theta,dtheta,T,dT,u,alpha,Vl,delta_u,para)
% The 'Linearization_B' function is to compute the approximate system matrix 
% whose form is linear.
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
f = Aeq\Beq;

u = u+delta_u;
Beq_delta = zeros(3,1);
Beq_delta(1,1) = para.m*Vl*dT*cos(theta)-2*para.m*dtheta*Vl*dT*sin(theta)-...
    para.m*dtheta^2*Vl*T*cos(theta)-(para.rho*para.L^2/2+para.m*para.L)*dtheta^2*sin(theta)-...
    2*para.rho*dtheta*alpha(1)*dT*sin(theta)-para.rho*alpha(1)*dtheta^2*T*cos(theta)-u;
Beq_delta(2,1) = 2*para.m*dtheta*Vl^2*T*dT-para.m*para.g*Vl*T*cos(theta)+...
    2*dtheta*para.rho*alpha(2)*T*dT-para.rho*para.g*T*cos(theta)-...
    (para.m*para.g*para.L+para.rho*para.g*para.L^2/2)*sin(theta);
Beq_delta(3,1) = -para.m*Vl*T*dtheta^2-para.m*para.g*sin(theta)-alpha(4)*T;

f_delta = Aeq\Beq_delta;

dfdu = (f_delta-f)/delta_u;
B_approx = [0;dfdu(1);0;dfdu(2);0;dfdu(3)];