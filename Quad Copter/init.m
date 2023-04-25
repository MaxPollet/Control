%% Syms initialization

% states
syms x y z vx vy vz phi theta psi wx wy wz ;

% inputs
syms v21 v22 v23 v24

% parameters
syms m L k b g kd Ixx Iyy Izz cm;


%% for values

% m = 0.5;
% L = 0.25;
% k = 3*10^(-6);
% b = 10^(-7);
% g = 9.81;
% kd = 0.25;
% Ixx = 5*10^(-3);
% Iyy = 5*10^(-3);
% Izz = 10^(-2);
% cm = 10^(4);



%% Equilibrium point at y* = 0
% It is known through filling in xdot = 0 at y* = 0 that the equilibrium
% point is the following: 
%   x* = zeros(12,1)
%   u* = g*m/(4k*cm)*ones(4,1)

x_0 = 0;
y_0 = 0;
z_0 = 0;
vx_0 = 0;
vy_0 = 0;
vz_0 = 0;
psi_0 = 0;
phi_0 = 0;
theta_0 = 0;
wx_0 = 0;
wy_0 = 0;
wz_0 = 0;

v21_0 = g*m/(4*k*cm);
v22_0 = g*m/(4*k*cm);
v23_0 = g*m/(4*k*cm);
v24_0 = g*m/(4*k*cm);
%% Nonlinear state equations of the system

f1 = vx;
f2 = vy;
f3 = vz;
f4_a = -(kd/m)*vx; 
f4_b = (k*cm/m)*(sin(psi)*sin(phi)+cos(phi)*cos(psi)*sin(theta))*(v21+v22+v23+v24);

f4 = f4_a + f4_b;

f5_a = -(kd/m)*vy;
f5_b = (k*cm/m)*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi))*(v21+v22+v23+v24);

f5 = f5_a + f5_b;

f6_a = -(kd/m)*vz - g;
f6_b = (k*cm/m)*(cos(theta)*cos(phi))*(v21+v22+v23+v24);

f6 = f6_a + f6_b;

f7 = wx + wy*(sin(phi)*tan(theta)) + wz*(cos(phi)*tan(theta));
f8 = wy*cos(phi) - wz*sin(phi);
f9 = sin(phi)/cos(theta)*wy + cos(phi)/cos(theta)*wz;
f10_a =  -((Iyy- Izz)/Ixx)*wy*wz;
f10_b = (L*k*cm/Ixx)*(v21- v23);

f10 = f10_a + f10_b;

f11_a = -((Izz- Ixx)/Iyy)*wx*wz;
f11_b = (L*k*cm/Iyy)*(v22- v24);

f11 = f11_a + f11_b;

f12_a =  -((Ixx - Iyy)/Izz)*wx*wy;
f12_b = (b*cm/Izz)*(v21- v22+ v23 - v24);

f12 = f12_a + f12_b;

%% Linearization of the state space equations

% matrix A
J_a = jacobian([f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12],[x y z vx vy vz phi theta psi wx wy wz]);

% matrix B
J_b = jacobian([f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12],[v21 v22 v23 v24]);

%% Evaluating the jacobian in x* and u*

Js_a = subs(J_a, [x y z vx vy vz phi theta psi wx wy wz v21 v22 v23 v24], ...
    [x_0 y_0 z_0 vx_0 vy_0 vz_0 phi_0 theta_0 psi_0 wx_0 wy_0 wz_0 v21_0 v22_0 v23_0 v24_0]); 

Js_b = subs(J_b, [x y z vx vy vz theta psi phi wx wy wz v21 v22 v23 v24], ...
    [x_0 y_0 z_0 vx_0 vy_0 vz_0 phi_0 theta_0 psi_0 wx_0 wy_0 wz_0 v21_0 v22_0 v23_0 v24_0]);



% Filling in the parameters

A = double(subs(Js_a, [m L k b g kd Ixx Iyy Izz cm ], ...
    [0.5 0.25 3*10^(-6) 10^(-7) 9.81 0.25 5*10^(-3) 5*10^(-3) 10^(-2) 10^(4) ])); 
B = double(subs(Js_b, [m L k b g kd Ixx Iyy Izz cm], ...
    [0.5 0.25 3*10^(-6) 10^(-7) 9.81 0.25 5*10^(-3) 5*10^(-3) 10^(-2) 10^(4)])); 
C1 = [eye(3), zeros(3), zeros(3), zeros(3); zeros(3), zeros(3), eye(3), zeros(3)];
C = eye(12);
D = zeros(12,4);

D1 = zeros(6,4);

System = ss(A,B,C,D);

vstar = double(subs(g*m/(4*k*cm), [g m k cm], [9.81 0.5 3e-6 1e4]));
ystar = 0;

fprintf("The poles of the linearized system:\n")
pole(System)
%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization
%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 0.05;

fprintf("Zero hold discretization\n")
Zero_hold = c2d(System,Ts);
fprintf("The poles:\n")
pole(Zero_hold)
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(Zero_hold)))
fprintf("The rank of the observability matrix: %d\n", rank(obsv(Zero_hold)))

fprintf("Bilinear discretization\n")
Tustin = c2d(System,Ts, 'tustin');
fprintf("The poles:\n")
pole(Tustin)
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(Tustin)))
fprintf("The rank of the observability matrix: %d\n", rank(obsv(Tustin)))

fprintf("Forward Euler discretization\n")
euler = c2d(System,Ts, 'forward');
fprintf("The poles:\n")
pole(euler)
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(euler)))
fprintf("The rank of the observability matrix: %d\n", rank(obsv(euler)))


[a,b,c,d] = ssdata(Tustin); 

%%%%%%%%%%%%%%%%%%%%%%%%
%% LQG control
%%%%%%%%%%%%%%%%%%%%%%%%



m = 4;
n = 12;


% %%%%%%%%%%%%%%
% %INPUT > OUTPUT;
% %%%%%%%%%%%%%%%

Q = eye(12)*10;
Q(1,1) = 50;
Q(2,2) = 50;
Q(3,3) = 400;

R = eye(4)*0.005;

[K,S,P] = dlqr(a,b,Q,R);



L = [a - eye(12), b];
O = [c,d];
I =  [L; O];


Z = [zeros(n,12); eye(12)];

N = I\Z;


Nx = N(1:n,:);
Nu = N(end-(m-1):end,:);










