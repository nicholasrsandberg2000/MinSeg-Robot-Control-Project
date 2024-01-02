% Nicholas Sandberg
% 015676752
% ME190 Section 4
% December 5, 2023
% Minseg Control Project

%% Part 1: Modeling the MinSeg Robot
clear all; close all; clc;
% Parameters (with units)
mp = 0.285; % kg
mw = 0.025; % kg
Jp = 0.0010; % kg*m^2
Jw = 0.0013; % kg*m^2
rw = 0.022; % m
lp = 0.1; % m
c = 0.0001; % N*m*s/rad
Kt = 0.3; % N*m/A
Kb = 0.5; % V*s/rad
R = 5; % ohms
g = 9.81; % m/s^2

% Linearized equation in matrix form coefficients
M = [((Jw/rw^2)+mw+mp),(mp*lp);(mp*lp),(mp*(lp^2)+Jp)];
D = [(c/rw^2+(Kb*Kt)/(R*rw^2)),((-Kb*Kt)/(R*rw)-(c/rw));((-Kb*Kt)/(R*rw)-(c/rw)),(((Kb*Kt)/R)+c)];
K = [0,0;0,(-mp*g*lp)];
F = [(Kt/(R*rw));(-Kt/R)];

% Convert matrix form equation to state-space model
a1 = zeros(2); 
a2 = eye(2);
a3 = -inv(M)*K;
a4 = -inv(M)*D;

b1 = zeros(2,1);
b2 = inv(M)*F;

A = [a1,a2;a3,a4];
B = [b1;b2];

% Verify eigenvalues of A and values of matrix B. eig(A) = [0;6.9;-6.0;-38] B = [0;0;1.14;-24.01]
a = eig(A)
b = B

% Simulate state-space system for small initial angle for zero input
% voltage
theta1 = pi/60;
q0 = [0;theta1;0;0];

% Generate output equation
C1 = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
D = [0;0;0;0];

% Run state space simulation and plot results
Sys = ss(A,B,C1,D);
[y,t,q] = initial(Sys,q0,0.1);

% Wheel Position vs. Time
figure(1);
plot(t,q(:,1),'r','linewidth',2)
title('Open-loop System')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(2); 
plot(t,q(:,2),'g','linewidth',2)
title('Open-loop System')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(3); 
plot(t,q(:,3),'b','linewidth',2)
title('Open-loop System')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(4); 
plot(t,q(:,4),'m','linewidth',2)
title('Open-loop System')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

%% Part 2: Designing a Full-State Feedback Controller
% Desired location of poles
P_des = [-6 -6 -6 -6];

% Gain calculation using Ackerman method
K_acker = acker(A,B,P_des)

% A matrix based on calculated gain values and simulate state space system
A_cl = (A-B*K_acker);
Sys_CL = ss(A_cl,B,C1,D);

Eig_acker = eig(A_cl)

% Run state space simulation and plot results
[y,t,x] = initial(Sys_CL,q0);

% Wheel Position vs. Time
figure(5); plot(t,x(:,1),'r','Linewidth',2)
title('Closed-loop System: Pole Placement with Ackerman method ')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(6); plot(t,x(:,2),'g','Linewidth',2)
title('Closed-loop System: Pole Placement with Ackerman method ')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(7); plot(t,x(:,3),'b','Linewidth',2)
title('Closed-loop System: Pole Placement with Ackerman method ')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(8); plot(t,x(:,4),'m','Linewidth',2)
title('Closed-loop System: Pole Placement with Ackerman method ')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
theta2 = (10*pi/180);
q0 = [0;theta2;0;0];
u = -K_acker*x';
figure(9); plot(t,u,'k','Linewidth',2)
title('Closed-loop System: Pole Placement with Ackerman method ')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

%% Part 3: Designing a Linear Quadratic Regulator (LQR)
% Define weight matrices
Q = diag([1,100,1,100]);
R = 1;

% Obtain gain value using LQR function
K_LQR = lqr(A,B,Q,R)

% A matrix based on calculated gain values and simulate state space system
A_cl = A-B*K_LQR;
Sys_CL_LQR = ss(A_cl,B,C1,D);

Eig_LQR = eig(A_cl)

% Run state space simulation and plot results
[y,t,x] = initial(Sys_CL_LQR,q0);

% Wheel Position vs. Time
figure(10); plot(t,x(:,1),'r','Linewidth',2)
title('Closed-loop control: LQR')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(11); plot(t,x(:,2),'g','Linewidth',2)
title('Closed-loop control: LQR')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(12); plot(t,x(:,3),'b','Linewidth',2)
title('Closed-loop control: LQR')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(13); plot(t,x(:,4),'m','Linewidth',2)
title('Closed-loop control: LQR')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
theta3 = (10*pi/180);
q0 = [0;theta3;0;0];
u = -K_LQR*x';
figure(14); plot(t,u,'k','Linewidth',2)
title('Closed-loop control: LQR')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

%% Part 4.8 Gain Calculation
close all; clc;
% Define weight matrices
Q = diag([100,100,1,100]);
R = 1;

% Obtain gain value using LQR function
K_LQR = lqr(A,B,Q,R)

%% Part 4.8 Plots
% Wheel Position vs. Time
figure(1); plot(t,x,'r','Linewidth',2)
title('Feedback Controller: Identical weights for x and theta')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(2); plot(t,Theta,'g','Linewidth',2)
title('Feedback Controller: Identical weights for x and theta')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(3); plot(t,x_dot,'b','Linewidth',2)
title('Feedback Controller: Identical weights for x and theta')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(4); plot(t,Theta_dot,'m','Linewidth',2)
title('Feedback Controller: Identical weights for x and theta')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
figure(5); plot(t,V,'k','Linewidth',2)
title('Feedback Controller: Identical weights for x and theta')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

%% Part 4.9 Gain Calculation
clc;
% Define weight matrices
Q = diag([1,100,1,100]);
R = 1;

% Obtain gain value using LQR function
K_LQR = lqr(A,B,Q,R)

A_cl = A-B*K_LQR;
Sys_CL_LQR = ss(A_cl,B,C1,D);

Eig_LQR = eig(A_cl)

%% Part 4.9 Plots
close all;
% Wheel Position vs. Time
figure(1); plot(t,x,'r','Linewidth',2)
title('Feedback Controller Data')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(2); plot(t,Theta,'g','Linewidth',2)
title('Feedback Controller Data')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(3); plot(t,x_dot,'b','Linewidth',2)
title('Feedback Controller Data')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(4); plot(t,Theta_dot,'m','Linewidth',2)
title('Feedback Controller Data')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
figure(5); plot(t,V,'k','Linewidth',2)
title('Feedback Controller Data')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

%% Part 4.10 Accelerometer Only
close all;
% Wheel Position vs. Time
figure(1); plot(t,x,'r','Linewidth',2)
title('Feedback Controller Data: Accelerometer Only')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(2); plot(t,Theta,'g','Linewidth',2)
title('Feedback Controller Data: Accelerometer Only')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(3); plot(t,x_dot,'b','Linewidth',2)
title('Feedback Controller Data: Accelerometer Only')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(4); plot(t,Theta_dot,'m','Linewidth',2)
title('Feedback Controller Data: Accelerometer Only')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
figure(5); plot(t,V,'k','Linewidth',2)
title('Feedback Controller Data: Accelerometer Only')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

%% Part 4.10 Gyroscope Only
close all;
% Wheel Position vs. Time
figure(1); plot(t,x,'r','Linewidth',2)
title('Feedback Controller Data: Gyroscope Only')
xlabel('Time (s)')
ylabel('Wheel Position (m)')
legend('x1 (m)')

% Pendulum Angle vs. Time
figure(2); plot(t,Theta,'g','Linewidth',2)
title('Feedback Controller Data: Gyroscope Only')
xlabel('Time (s)')
ylabel('Pendulum Angle (rad)')
legend('Theta (rad)')

% Wheel Velocity vs. Time
figure(3); plot(t,x_dot,'b','Linewidth',2)
title('Feedback Controller Data: Gyroscope Only')
xlabel('Time (s)')
ylabel('Wheel Velocity (m/s)')
legend('x2 (m/s)')

% Pendulum Angle Vector vs. Time
figure(4); plot(t,Theta_dot,'m','Linewidth',2)
title('Feedback Controller Data: Gyroscope Only')
xlabel('Time (s)')
ylabel('Pendulum Angle Vector (rad/s)')
legend('Theta dot (rad/s)')

% Plot step response with initial condition of 10 degrees
figure(5); plot(t,V,'k','Linewidth',2)
title('Feedback Controller Data: Gyroscope Only')
xlabel('Time (s)')
ylabel('Input Voltage (V)')
legend('Input voltage (V)')

