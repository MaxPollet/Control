syms vx vy vz theta psi phi wx wy wz v21 v22 v23 v24 x y z;


%% Voor puur in symbolen
syms m L k b g kd Ixx Iyy Izz cm;


%% voor een uitgewerkte versie

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

%% equi point values
psi_0 = 0;
phi_0 = 0;
theta_0 = 0;
x_0 = 0;
y_0 = 0;
z_0 = 0;

f1 = vx;
f2 = vy;
f3 = vz;
f4_a = -(kd/m)*vx; 
f4_b = (kd*cm/m)*(sin(psi)*sin(phi)+cos(phi)*cos(psi)*sin(theta))*(v21+v22+v23+v24);
f5_a = -(kd/m)*vy;
f5_b = (kd*cm/m)*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(theta))*(v21+v22+v23+v24);
f6_a = -(kd/m)*vz -g;
f6_b = (kd*cm/m)*(cos(theta)*cos(phi))*(v21+v22+v23+v24);
f7 = wx + wy*(sin(phi)*tan(theta)) + wz*(cos(phi)*tan(theta));
f8 = wy*cos(phi) - wz*sin(phi);
f9 = sin(phi)/cos(theta)*wy + cos(phi)/cos(theta)*wz;
f10_a =  -((Iyy- Izz)/Ixx)*wy*wz;
f10_b = (L*k*cm/Ixx)*(v21- v23);
f11_a = -((Izz- Ixx)/Iyy)*wx*wz;
f11_b = (L*k*cm/Iyy)*(v22- v24);
f12_a =  -((Ixx - Iyy)/Izz)*wx*wy;
f12_b = (L*k*cm/Izz)*(v21- v21+ v24 - v24);

%% jacobiaan om de a matrix te maken
J_a = jacobian([f1 f2 f3 f4_a f5_a f6_a f7 f8 f9 f10 f11 f12],[vx vy vz theta phi psi x y z wx wy wz]);

%% jacobiaan om de B matrix te maken
J_b = jacobian([0 0 0 f4_b f5_b f6_b 0 0 0 f10_b f11_b f12_b],[v21 v22 v23 v24]);

%% De equilibrium waardes in de Jacobiaan steken
%maar omdat er geen wx wy wz equib waardes zijn is het nog steeds
%nonlineair???
Js_a = subs(J_a, [x y z theta psi phi], [0 0 0 0 0 0]); 
%A = double(Js_a);
Js_b = subs(J_b, [x y z theta psi phi], [0 0 0 0 0 0]); 
%B = double(Js_b);


System = ss(Js_a,Js_b);