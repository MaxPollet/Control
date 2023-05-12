clc
clear all
close all
%% Parameters init
K_stiff = 1.60856;  % Spring stiffness
K_g     = 70;       % Motor gear ratio
K_m     = 0.00767;  % Motor torque constant
K_b     = 0.00767;  % Motor back EMF constant
J_h     = 0.0021;   % Hub inertia
R_m     = 2.6;      % Motor armature resistance
J_l     = 0.0059;   % Arm inertia
d       = 0.0318;   % Body anchor point (y)
R       = 0.076;    % Arm anchor point
r       = 0.0318;   % Body anchor point (x)
F_r     = 1.33;     % Spring restoring force
L       = 0.0318;   % Spring length at rest

%% State space description
A = [0,  0,                             1,                           0; 
     0,  0,                             0,                           1; 
     0,  K_stiff/J_h,                  -((K_g)^2*K_m*K_b)/(J_h*R_m), 0;
     0, -((J_h+J_l)*K_stiff)/(J_h*J_l), ((K_g)^2*K_m*K_b)/(J_h*R_m), 0];

B = [ 0; 
      0; 
      (K_m*K_g)/(R_m*J_h);
     -(K_m*K_g)/(R_m*J_h)];

C = [1 0 0 0;
     0 1 0 0];
D = zeros(2,1);

System = ss(A, B, C, D);
TF = tf(System);
fprintf("The eigenvalues of the matrix A:")
eig(A)
fprintf("The poles of the system are:")
pole(System)
disp("transmission zeros of the system")
disp(tzero(A, B, C, D))
fprintf("The rank of the controllability matrix is: %d\n", rank(ctrb(A,B)));
fprintf("The rank of the observability matrix is: %d\n", rank(obsv(A,C)));

figure
pzmap(System)
%% Controller

Q = {[350, 0,    0, 0; 
       0,   1500, 0, 0; 
       0,   0,    3, 0;
       0,   0,    0, 0.5]; % initiele controller
       [400, 0,    0, 0; 
       0,   1500, 0, 0; 
       0,   0,    3, 0;
       0,   0,    0, 0.5]; % dit lijkt me de optimale
       [350, 0,    0, 0; 
       0,   1000, 0, 0; 
       0,   0,    3, 0;
       0,   0,    0, 0.5];
       [350, 0,    0, 0; 
       0,   1500, 0, 0; 
       0,   0,    50, 0;
       0,   0,    0, 0.5]; % voorbeeld van te trage controller
       [500, 0,    0, 0; 
       0,   1500, 0, 0; 
       0,   0,    3, 0;
       0,   0,    0, 0.5]; % voorbeeld van te nerveuze controller
       [350, 0,    0, 0; 
       0,   1500, 0, 0; 
       0,   0,    3, 0;
       0,   0,    0, 0.5]};

R = {10 10 10 10 0.05 15};

for i=1:length(Q)
    Q_c = Q{i};
    R_c = R{i};
    
    [K, ~, P] = lqr(System, Q_c, R_c);
    
    %out = sim("controller_check_2021.slx");
    
%     figure
%     subplot(2, 1, 1)
%     plot(out.tout, out.theta.Data, 'Color', "#0072BD")
%     hold on
%     plot(out.tout, out.theta_ref.Data, 'Color', "#D95319")
%     grid on
%     xlabel("$$t [s]$$", 'Interpreter','latex')
%     ylabel("$$\theta [rad]$$", 'Interpreter', 'latex')
%     legend('\theta', '\theta_{ref}')
%     subplot(2,1,2)
%     plot(out.tout, out.alpha.Data, 'Color', "#0072BD")
%     grid on
%     xlabel("$$t [s]$$", 'Interpreter','latex')
%     ylabel("$$\alpha [rad]$$", 'Interpreter', 'latex')
% 
%     figure
%     subplot(2, 1, 1)
%     plot(out.tout, out.theta_dot.Data, 'Color', "#0072BD")
%     grid on
%     xlabel("$$t [s]$$", 'Interpreter','latex')
%     ylabel("$$\dot\theta \left[\frac{1}{s}\right]$$", 'Interpreter', 'latex')
%     subplot(2,1,2)
%     plot(out.tout, out.alpha_dot.Data, 'Color', "#0072BD")
%     grid on
%     xlabel("$$t [s]$$", 'Interpreter','latex')
%     ylabel("$$\dot\alpha \left[\frac{1}{s}\right]$$", 'Interpreter', 'latex')
% 
%     figure
%     plot(out.tout, out.control_actions.Data, 'Color', "#0072BD")
%     grid on
%     xlabel("$$t [s]$$", 'Interpreter','latex')
%     ylabel("$$V [Volt]$$", 'Interpreter', 'latex')
% 
%     % werkt niet voor een reden
% %     bode(close_sys)
% 
%     disp("-------------------------")
%     disp(i)
%     disp("The closed loop eigevalues: ")
%     disp(P)
    
end

Q_c = Q{2};
R_c = R{2};

[K,S,P] = lqr(System,Q_c,R_c);

fprintf("The closed-loop eigenvalues are:")
display(P)

controller = tf(ss(A-B*K, zeros(4,1), C, zeros(2,1)));
% 
% figure
% bode(controller)

%%
%%%%%%%%%%%%%%%%
% Implementation
%%%%%%%%%%%%%%%%

w_c = 2*2*pi;
Ts = 1/200;

filter = tf([w_c*Ts,0], [1+w_c*Ts, -1]);


