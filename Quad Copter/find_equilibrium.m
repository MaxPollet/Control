function [xstar, ustar] = find_equilibrium(ystar, param)
% This function finds the value of x* and u* given y* and the parameters
    %% initilization
    xstar = zeros(12,1);
    ustar = zeros(4,1);

    x     = ystar(1);
    y     = ystar(2);
    z     = ystar(3);
    phi   = ystar(4);
    theta = ystar(5);
    psi   = ystar(6);

    m   = param(1);
    L   = param(2);
    k   = param(3);
    b   = param(4);
    g   = param(5);
    kd  = param(6);
    Ixx = param(7);
    Iyy = param(8);
    Izz = param(9);
    cm  = param(10);
    %% filling in the known states

    xstar(1:3) = [x; y; z];
    xstar(7:9) = [phi; theta; psi];
    
    %% Solving for unknown states and inputs

    F = @(x) [k*cm/m*(sin(psi)*sin(phi) + cos(psi)*cos(phi)*sin(theta))*(x(4) + x(5) + x(6) + x(7));
              k*cm/m*(sin(psi)*cos(phi)*sin(theta) - cos(psi)*sin(phi))*(x(4) + x(5) + x(6) + x(7));
              -g + k*cm/m*(cos(theta)*cos(phi))*(x(4) + x(5) + x(6) + x(7));
              x(1) + x(2)*(sin(phi)*tan(theta)) + x(3)*(cos(phi)*tan(theta));
              x(2)*cos(phi) - x(3)*sin(phi);
              sin(phi)/cos(theta)*x(2) + cos(phi)/cos(theta)*x(3);
              L*k*cm/Ixx*(x(4) - x(5)) - (Iyy - Izz)/Ixx*x(2)*x(3);
              L*k*cm/Iyy*(x(5) - x(7)) - (Izz - Ixx)/Iyy*x(1)*x(2);
              L*k*cm/Izz*(x(4) - x(5) + x(6) - x(7)) - (Ixx - Iyy)/Izz*x(1)*x(2)];
    x0 = zeros(7,1);
    x = fsolve(F, x0);

    xstar(10:12) = x(1:3);
    ustar = x(4:7);
end