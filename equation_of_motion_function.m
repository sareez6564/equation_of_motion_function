function [T1,T2,T3] = equation_of_motion_function(q, qdot, qddot, A)
    
% Define constants
m = [0; 1.17; 0.83];
L = [0; 0.7; 0.5];
r = [0; 0.02; 0.02];
p = [0; 1330.0806; 1320.98603];

% Distance from each frame to the associated body's com (respect to b-f-f)
Lc1 = [0; -L(1)/2; 0];
Lc2 = [-L(2)/2; 0; 0];
Lc3 = [-L(3)/2; 0; 0];

% Define Inertia Matrices with respect to the bode-fixed-frame
I1=[p(1)*((2/3)*pi*r(1)*(L(1)^3)+(1/3)*pi*(r(1)^3)*L(1)) 0 0;
    0 p(1)*((2/3)*pi*L(1)*(r(1)^3)) 0;
    0 0 p(1)*((1/3)*pi*L(1)*(r(1)^3)+(2/3)*pi*(r(1))*(L(1)^3))];

I2 = [p(2)*((2/3)*pi*(r(2)^3)*L(2)), 0, 0;
          0, p(2)*((2/3)*pi*r(2)*(L(2)^3) + (1/3)*pi*(r(2)^3)*L(2)), 0;
          0, 0, p(2)*((2/3)*pi*r(2)*(L(2)^3) +(1/3)*pi*(r(2)^3)*L(2))];

I3 = [p(3)*((2/3)*pi*(r(3)^3)*L(3)), 0, 0;
          0, p(3)*((2/3)*pi*r(3)*(L(3)^3) + (1/3)*pi*(r(3)^3)*L(3)), 0;
          0, 0, p(3)*((2/3)*pi*r(3)*(L(3)^3) +(1/3)*pi*(r(3)^3)*L(3))];

% Rotation matrix from frame 0 to 1
R1 = [cos(q(1)) 0 sin(q(1));
          sin(q(1)) 0 -cos(q(1));
          0 1 0];

% Rotation matrix from frame 1 to 2
R2 = [cos(q(2)) -sin(q(2)) 0;
          sin(q(2)) cos(q(2)) 0;
          0 0 1];

% Rotation matrix from frame 2 to 3
R3 = [cos(q(3)) -sin(q(3)) 0;
          sin(q(3)) cos(q(3)) 0;
          0 0 1];
          
% Forward Recursion
w0 = [0;0;0];     % Base has zero angular velocity
w0dot = [0;0;0];  % Base has zero angular acceleration
p0ddot = [A;0;9.81]; % Linear accelearion of the base
z0 = [0;0;1];     % Rotation about the base z-axis

% Angular velocities of each link with respect to their associated b-f-f
w1 = R1' * (w0 + qdot(1)*z0);
w2 = R2' * (w1 + qdot(2)*z0);
w3 = R3' * (w2 + qdot(3)*z0);

% Angular accelerations of each link with respect to their associated b-f-f
w1dot = R1' * (w0dot + qddot(1)*z0 + qdot(1)*cross(w0,z0));
w2dot = R2' * (w1dot + qddot(2)*z0 + qdot(2)*cross(w1,z0));
w3dot = R3' * (w2dot + qddot(3)*z0 + qdot(3)*cross(w2,z0));

% Distance from frame 0 to frame 1 with respect to frame 1(respect to F1)
r0to1 = [0;L(1);0];
% Distance from frame 1 to frame 2 with respect to frame 2(respect to F2)
r1to2= [L(2);0;0];
% Distance from frame 2 to frame 3 with respect to frame 3(respect to F3)
r2to3= [L(3);0;0];
 
% Linear accelerations of each link with respect to their associated b-f-f
p1ddot = R1'*p0ddot + cross(w1dot,r0to1) + cross(w1,cross(w1,r0to1));
p2ddot = R2'*p1ddot + cross(w2dot,r1to2) + cross(w2,cross(w2,r1to2));
p3ddot = R3'*p2ddot + cross(w3dot,r2to3) + cross(w3,cross(w3,r2to3));

% Linear accelerations of each link's com with respect to their associated b-f-f
p1ddotc= p1ddot + cross(w1dot,Lc1) + cross(w1,cross(w1,Lc1));
p2ddotc= p2ddot + cross(w2dot,Lc2) + cross(w2,cross(w2,Lc2));
p3ddotc= p3ddot + cross(w3dot,Lc3) + cross(w3,cross(w3,Lc3));

% Backward Recursion
% Forces transmitted onto each body with respect to their associated b-f-f
f3 = m(3)*p3ddotc;
f2 = R3*f3 + m(2)*p2ddotc;
f1 = R2*f2 + m(1)*p1ddotc;

% Moments transmitted onto each body with respect to their associated b-f-f
tau3 = -cross(f3, (r2to3 + Lc3)) + I3 * w3dot + cross(w3, I3 * w3);
tau2 = -cross(f2, (r1to2 + Lc2)) + R3 * tau3 + R3 * cross(f3, Lc2) + I2 * w2dot + cross(w2, I2 * w2);
tau1 = -cross(f1, (r0to1 + Lc1)) + R2 * tau2 + R2 * cross(f2, Lc1) + I1 * w1dot + cross(w1, I1 * w1);

% Compute generalized torques
T1 = tau1' * R1' * z0;
T2 = tau2' * R2' * z0;
T3 = tau3' * R3' * z0;

end
