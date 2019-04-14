%% Simple Descent Simulation: Design Spring 2019

% Author: Henry Macanas
% Version: 1
% Date: 3/24/19

% housekeeping ````````````````````````````````````````````````````````````

clc, close all, clear all

% constants ```````````````````````````````````````````````````````````````

const = struct;                          % structure of constant
const.mu = 4.9048695e12;                 % m3/s2, grav param moon
const.m = 2000;                          % kg, initilal mass of spacecraft
const.len = 3;                           % 3m, length of box
const.I = diag(1/6*const.m*[1; 1; 1]);   % kg-m2, inertia of spacecraft
const.r_moon = 1737.1 * 1000;            % m, radius of the moon
const.u_body = [0; 0; -1];               % main engine thrust vect in body

% initial conditions `````````````````````````````````````````````````````` 

r_0 = [0; 0; 1000 + const.r_moon];       % m, position vector tbd
v_0 = [0; 0; 0];                         % m/s, velocity vector tbd
eul_0 = [deg2rad(45); 0; 0];             % rads, euler vector tbd
w_0 = [0; 0; 0];                         % rad/s, euler vector tbd
state = [r_0; v_0; eul_0; w_0];          % states of system tbd

% call ode ````````````````````````````````````````````````````````````````

tspan = [0 60];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[tnew, statenew] = ode45(@sc_dynamics,tspan,state,options,const);

% plots ```````````````````````````````````````````````````````````````````

% descent profile 3D
figure
plot3(statenew(:, 1), statenew(:, 2), statenew(:, 3));
title("PDV Descent Profile")
xlabel("x - pos (m)");
ylabel("y - pos (m)");
zlabel("z - pos (m)");
grid on

% descent profile yz
figure
plot(statenew(:, 2), statenew(:, 3));
hold on
plot(statenew(end, 2) ,statenew(end, 3),'r*')
plot(statenew(1, 2), statenew(1, 3),'g*')
title("PDV Descent Profile y-z")
xlabel("y - pos (m)");
ylabel("z - pos (m)");
grid on

% descent dynamics ` ``````````````````````````````````````````````````````

function [y] = sc_dynamics(t, state, const)

% define state parameters `````````````````````````````````````````````````
r = state(1:3);                                 % m, pos vect in inertial
v = state(4:6);                                 % m/s, velocity vector
eul_b = state(7:9);                             % rads, euler angles
w_b= state(10:12);                              % rad/s, ang vel body

% transformations `````````````````````````````````````````````````````````
C = cx(eul_b(1)) * cy(eul_b(2)) * cz(eul_b(3)); % body to inertial 
r_b = C * r;                                    % m, pos vect in body
u_i = C' * const.u_body;                        % main eng thrust vect iner

% disturbance torques `````````````````````````````````````````````````````

r_mag = norm(r);                                % m, magnitude of postion 
T_g = 3*const.mu/r_mag^5*cross_matrix(r_b)*const.I*r_b;

% eoms ````````````````````````````````````````````````````````````````````

T = 1000;                                       % main engine thrust
dv = -const.mu*r/r_mag^3 + T/const.m*u_i;       % m/s2, current acc vector
dw = const.I\(T_g - cross(w_b,const.I*w_b));
y = [v; dv; w_b; dw];

end

