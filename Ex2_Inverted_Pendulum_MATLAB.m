clear global
close all
clear
clc

global M m L g im
M = 0.5;      % mass of the cart (kg)  
m = 0.5;      % mass of the pendulum bob (kg)
L = 1;        % length of the pendulum rod (m)
g = 9.81;
im = 0.5;      % inertia  

a21 = (m*g*(L/2))/(im+m*(L/2)^2-((m^2*(L/2)^2)/(m+M)));
a41 = ((im+m*(L/2)^2)*g)*(im+m*(L/2)^2-(m^2*(L/2)^2)/(M+m));
b2 = ((m*(L/2))/(M+m))/(im+m*(L/2)^2-(m^2*(L/2)^2)/(M+m));
b4 = ((im+m*(L/2)^2)/(M+m))/(im+m*(L/2)^2-(m^2*(L/2)^2)/(M+m));

A = [0 1 0 0;
    a21 0 0 0;
    0 0 0 1;
    a41 0 0 0];
B = [0;
    b2;
    0;
    b4];
ns = size(A,2);
ms = size(B,2);
%% Variable define
X=sdpvar(ns,ns,'symmetric');
F = sdpvar(ms,ns); 
LMI1= X*A'+F'*B'+A*X+B*F<=0;

LMI=LMI1+[X>=0];
ops=sdpsettings('solver','sedumi','showprogress',0,'verbose',0);
solvesdp(LMI,[],ops)
checkset(LMI)

K = double(F)*inv(double(X));

%% Linearization about x1=0, LQR,  Simulation
tf=10;  % final time
ti=0.001;  % runge kutta sample time 
tspan=0:ti:tf;
sample_size = size(tspan,2);

% x(:,1)=[0.1745;0;0;0]; 
degree = 37
x(:,1)=[degree*pi/180;0;0;0]; 
for i=1:sample_size-1

    U = K*x(:,i);
    x(:,i+1) = rk6(x(:,i),U,ti);
    
end
figure(1);
plot(tspan,x(1,:)','r');
hold on
% plot(tspan,x(2,:)','b');
% hold on 
% plot(tspan,x(3,:)','g-');
% hold on 
% plot(tspan,x(4,:)','m-');
xlabel('Time (sec)');
ylabel('State');
legend('x(1)','x(2)','x(3)','x(4)')
grid on



function dx=Nplant(x,u)
global M m L g im
dx(1,1) = x(2);
dx(2,1) = (-m^2*(L/2)^2*sin(x(1))*cos(x(1))*x(2)^2+(M+m)*m*g*(L/2)*sin(x(1))+m*(L/2)*cos(x(1))*u)/((M+m)*(im+m*(L/2)^2)-m^2*(L/2)^2*cos(x(1))^2);
dx(3,1) = x(4);
dx(4,1) = (m^2*(L/2)^2*sin(x(1))*cos(x(1))*g-m*(im+m*(L/2)^2)*(L/2)*sin(x(1))*x(2)^2+(im+m*(L/2)^2)*u)/((M+m)*(im+m*(L/2)^2)-m^2*(L/2)^2*cos(x(1))^2);
end
function dx=rk6(x,u,T)
k1=Nplant(x,u)*T;
k2=Nplant(x+k1*0.5,u)*T;
k3=Nplant(x+k2*0.5,u)*T;
k4=Nplant(x+k3,u)*T;
dx=x + ((k1+k4)/6+(k2+k3)/3);
end












