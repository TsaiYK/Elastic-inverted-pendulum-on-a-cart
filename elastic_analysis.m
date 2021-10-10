function [V,alpha] = elastic_analysis(theta0,seg_x,para)
% Reference:
% [1] Laura, P. A. A., J. L. Pombo, and E. A. Susemihl. "A note on the 
% vibrations of a clamped-free beam with a mass at the free end." 
% Journal of Sound and Vibration 37, no. 2 (1974): 161-168.
% Nguyen, Minh Hoang, Van Thuyen Ngo, Minh Tam Nguyen, Thi Thanh Hoang Le, 
% and V. D. H. Nguyen. "Designing linear feedback controller for elastic 
% inverted pendulum with tip mass." Robotica & Management 21, no. 2 (2016).

x_step = para.L/seg_x;
yl = 1.8751; % the root of the equation
q = para.rho*para.g*sin(theta0);
P = para.m*para.g*sin(theta0);
y_max = q*para.L^4/8/para.E/para.I+P*para.L^3/3/para.E/para.I;

R_V = ( cos(yl) + cosh(yl) )/( sin(yl) + sinh(yl) );
V_bar = cosh(yl*para.s/para.L) - R_V*sinh(yl*para.s/para.L) -...
    cos(yl*para.s/para.L) + R_V*sin(yl*para.s/para.L);
C1 = y_max/(cosh(yl) - R_V*sinh(yl) - cos(yl) + R_V*sin(yl));
% C1 = 0.5;

% Separation of variables
% v(s,t) = V(s)*T(t)
% space term, V(s)
V = C1*V_bar;


%% V(s) function
% coefficients alpha's
V_fun = @(x) (C1*cosh(yl*x/para.L) - R_V*C1*sinh(yl*x/para.L) -...
    C1*cos(yl*x/para.L) + R_V*C1*sin(yl*x/para.L));
alpha(1) = integral(V_fun,0,para.L);

V_sq_fun = @(x) (C1*cosh(yl*x/para.L) - R_V*C1*sinh(yl*x/para.L) -...
    C1*cos(yl*x/para.L) + R_V*C1*sin(yl*x/para.L)).^2;
alpha(2) = integral(V_sq_fun,0,para.L);

V_x_fun = @(x) x.*(C1*cosh(yl*x/para.L) - R_V*C1*sinh(yl*x/para.L) -...
    C1*cos(yl*x/para.L) + R_V*C1*sin(yl*x/para.L));
alpha(3) = integral(V_x_fun,0,para.L);

dVdx = (V(2:end)-V(1:end-1))/x_step;
ddVdx_sq = (dVdx(2:end)-dVdx(1:end-1))/x_step;
dddVdx_qd = (ddVdx_sq(2:end)-ddVdx_sq(1:end-1))/x_step;
alpha(4) = para.E*para.I*dddVdx_qd(end);
