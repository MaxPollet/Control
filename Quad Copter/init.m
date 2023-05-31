clear all
close all
clc
load("references_07.mat")
%% 
%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
disp("Linearizing the system...")
% states
syms x y z vx vy vz phi theta psi wx wy wz ;

% inputs
syms v21 v22 v23 v24

% parameters
syms m L k B_d g kd Ixx Iyy Izz cm;

% for values

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium point at y* = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear state equations of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
f12_b = (B_d*cm/Izz)*(v21- v22+ v23 - v24);

f12 = f12_a + f12_b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization of the state space equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matrix A
J_a = jacobian([f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12],[x y z vx vy vz phi theta psi wx wy wz]);

% matrix B
J_b = jacobian([f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12],[v21 v22 v23 v24]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluating the jacobian in x* and u*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Js_a = subs(J_a, [x y z vx vy vz phi theta psi wx wy wz v21 v22 v23 v24], ...
    [x_0 y_0 z_0 vx_0 vy_0 vz_0 phi_0 theta_0 psi_0 wx_0 wy_0 wz_0 v21_0 v22_0 v23_0 v24_0]); 

Js_b = subs(J_b, [x y z vx vy vz theta psi phi wx wy wz v21 v22 v23 v24], ...
    [x_0 y_0 z_0 vx_0 vy_0 vz_0 phi_0 theta_0 psi_0 wx_0 wy_0 wz_0 v21_0 v22_0 v23_0 v24_0]);


% Filling in the parameters

A_c = double(subs(Js_a, [m L k B_d g kd Ixx Iyy Izz cm ], ...
    [0.5 0.25 3*10^(-6) 10^(-7) 9.81 0.25 5*10^(-3) 5*10^(-3) 10^(-2) 10^(4) ])); 
B_c = double(subs(Js_b, [m L k B_d g kd Ixx Iyy Izz cm], ...
    [0.5 0.25 3*10^(-6) 10^(-7) 9.81 0.25 5*10^(-3) 5*10^(-3) 10^(-2) 10^(4)])); 
C_c = [eye(3), zeros(3), zeros(3), zeros(3); zeros(3), zeros(3), eye(3), zeros(3)];
D_c = zeros(6,4);

System = ss(A_c,B_c,C_c,D_c);

vstar = double(subs(g*m/(4*k*cm), [g m k cm], [9.81 0.5 3e-6 1e4]));
ystar = 0;

disp("The poles of the linearized system:")
disp(pole(System))

figure
pzmap(System)

clear x y z vx vy vz phi theta psi wx wy wz v21 v21_0 v22 v22_0 v23 v23_0 v24 v24_0 m L k B_d g kd Ixx Iyy Izz cm;
clear f1 f10 f10_a f10_b f11 f11_a f11_b f12 f12_a f12_b f2 f3 f4 f4_a f4_b f5 f5_a f5_b f6 f6_a f6_b f7 f8 f9;
clear J_a J_b Js_a Js_b 
%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Discretization
%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
fprintf("Beginning the discretization ...\n")
Ts = 0.05;

disp("Zero hold discretization")
Zero_hold = c2d(System,Ts);
disp("The poles:")
disp(pole(Zero_hold))
disp("The transmission zeros:")
disp(tzero(Zero_hold))
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(Zero_hold)))
fprintf("The rank of the observability matrix: %d\n\n", rank(obsv(Zero_hold)))

disp("Bilinear discretization")
Tustin = c2d(System,Ts, 'tustin');
disp("The poles:")
disp(pole(Tustin))
disp("The transmission zeros:")
disp(tzero(Tustin))
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(Tustin)))
fprintf("The rank of the observability matrix: %d\n\n", rank(obsv(Tustin)))

disp("Forward Euler discretization")
euler = c2d(System,Ts, 'forward');
disp("The poles:")
disp(pole(euler))
disp("The transmission zeros:")
disp(tzero(euler))
fprintf("The rank of the controlability matrix: %d\n", rank(ctrb(euler)))
fprintf("The rank of the observability matrix: %d\n\n", rank(obsv(euler)))

figure
pzmap(Tustin)

[A_d,B_d,C_d,D_d] = ssdata(Tustin); 
clear euler Zero_hold
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LQR control with full state feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
fprintf("Computing the gain for the LQR controller with full state feedback...\n")
m = 4;
n = 12;

%%%%%%%%%%%%%%
% without load
%%%%%%%%%%%%%%

I =  [A_d - eye(12), B_d; 
      C_d,           D_d];

Z = [zeros(n,6); eye(6)];

N = I\Z;

Nx = N(1:n,:);
Nu = N(end-(m-1):end,:);

Q = eye(12)*10;
Q(1,1) = 50;
Q(2,2) = 50;
Q(3,3) = 400;

R = eye(4)*0.1;

[K_fsf,S,P] = dlqr(A_d,B_d,Q,R);
disp("The gain of full state feedback with LQR")
disp(K_fsf)

%%%%%%%%%%%
% with load
%%%%%%%%%%%

Q = eye(12)*10;
Q(1,1) = 50;
Q(2,2) = 50;
Q(3,3) = 1000;

R = eye(4)*0.01;

[K_fsf_l,S,P] = dlqr(A_d,B_d,Q,R);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% LQR with integral action
%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
fprintf("Computing the gain for the LQR controller with integral control action...\n")

C_int = C_d(1:3,:);

NA = [eye(3),     C_d(1:3,:);
      zeros(n,3), A_d];

NB = [D_d(1:3,:);
      B_d];

fprintf("Rank of the controllability matrix of the augmented system: %d\n", ...
        rank(ctrb(NA, NB)))

Q_int = eye(15);
Q_int(4,4) = 50;
Q_int(5,5) = 50;
Q_int(6,6) = 350;

R_int = eye(4)*0.001;

[full_K, ~, P] = dlqr(NA, NB, Q_int, R_int);

Ki = full_K(:, 1:3);
Ks = full_K(:, 4:end);

disp("The gain of LQR with integral control")
disp("K_i:")
disp(Ki)
disp("K_s:")
disp(Ks)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state feedback via Pole placement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
disp("Pole placement with state feedback control and state estimation")
% Since there are 12 states do we need 12 poles, where 2 are dominant and
% 10 are non-dominant

% We want the settling time to be smaller than 7 seconds
% ts -> 2s
% damping ratio is trial and error
ts = 3;
dr = 0.7;

wn = 4.6/(dr*ts);
alpha = -dr*wn;
beta = wn*sqrt(1-dr^2);

%cl_poles_cont = [alpha - beta*1i, alpha + beta*1i, ones(1, 4)*alpha*8,... 
%                ones(1, 4)*alpha*9, ones(1, 2)*alpha*10]';

cl_poles_cont = [alpha - beta*1i, alpha + beta*1i, ones(1, 4)*alpha*3.1,... 
                ones(1, 4)*alpha*3.2, ones(1, 2)*alpha*3.3]';



cl_poles_disc = exp(cl_poles_cont.*Ts);

% computing the gain K
K_pp = place(A_d, B_d, cl_poles_disc);

% Now we make the estimator gain L
% First check if the system is observable
disp("The rank of the observability matrix:")
disp(rank(obsv(A_d, C_d)))

% The poles of the estimator should be 2-5 times bigger than the closed
% loop poles of the controller
cl_poles_est_cont = cl_poles_cont*2;

cl_poles_est_disc = exp(cl_poles_est_cont*Ts);

L_pp = place(A_d', C_d', cl_poles_est_disc)';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LQG with integral control
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("------------------------------------------------")
disp("LQG with integral control")

B1 = eye(12);
noise_sys = ss(A_d,[B_d, B1], C_d, [D_d,zeros(6,12)]);
A_n = A_d;
B_n = [B_d,B1];
C_n = C_d;
D_n = [D_d,zeros(6,12)];

%%%%% testing Q_noise
Q_noise = diag([1.5e-8, 1.5e-8, 1.5e-8, 6.5e-6, 6.5e-6, 6.5e-6, 6.5e-9, 6.5e-9, 6.5e-9, 2e-12, 2e-12, 2e-12]);
%%%% working Q_noise
%Q_noise = diag([1.5e-8, 1.5e-8, 1.5e-8, 6.5e-6, 6.5e-6, 6.5e-6, 6.5e-9, 6.5e-9, 6.5e-9, 2e-12, 2e-12, 2e-12]);


R_noise = [2.5e-5*eye(3),  zeros(3);
           zeros(3), 7.57e-5*eye(3)];

[K_est, S, P1] = dlqr(A_d', C_d', B1*Q_noise*B1', R_noise');
L_k = K_est';

[KEST,L_kalman,P] = kalman(noise_sys,Q_noise,R_noise);



