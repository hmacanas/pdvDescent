%% Simple Descent Simulation: Design Spring 2019

% Author: Henry Macanas
% Version: 1
% Date: 3/24/19

% housekeeping ````````````````````````````````````````````````````````````

clc, close all, clear all

% test desired profile ````````````````````````````````````````````````````

% initial conditions
r_0 = [0; -200; 500];        % m, position vector tbd
v_0 = [0; 20; -30];                         % m/s, velocity vector tbd

% target conditions
r_t = [0; 0; 5];
v_t = [0; 0; -1];
a_t = [0; 0; 0];

% compute travel time
t_go = getTGo(r_0(3), r_t(3), v_0(3), v_t(3));
t_cutoff = 1;
t_start = 0;
t_curr = t_start;
t_end = t_go;
step = .1;

% get trajectory
r = zeros(3, ceil(t_go/step));

for i = 1:ceil((t_go/ step))
    [C_0, C_1, C_2] = getTrajectoryCoeffs(r_0, r_t, v_0, v_t, a_t, t_go, t_cutoff);
    [a_prof, v_prof, r_prof] = getStateProfile(C_0, C_1, C_2, t_curr, v_0, r_0);
    r(:, i) = r_prof;
    t_curr = t_curr + step;
    t_go = t_go - step;
    
end 

figure 
plot(r(2, :), r(3, :))
hold on
plot(r(2, end), r(3, end),'r*')
plot(r(2, 1), r(3, 1),'g*')
title("PDV Descent Profile")
xlabel("x - pos (m)");
ylabel("y - pos (m)");
zlabel("z - pos (m)");
grid on
% functions ```````````````````````````````````````````````````````````````

function [t_go] = getTGo(r_0, r_t, v_0, v_t)

% inputs: 
    % r_0: vertical start postion (m)
    % r_t: target vertical end position (m)
    % v_0: vertical velocity start postion (m/s)
    % v_t: vertical velocity end position (m/s)
    % a_t: vertical target acceleration (m/s2)
    
 % output:
    % t_go: time till target reached (s)
    
    t_go = 3 * (r_t - r_0)/ (v_0 + 2 * v_t);
    
end

function [C_0, C_1, C_2] = getTrajectoryCoeffs(r_0, r_t, v_0, v_t, a_t, t_go, t_cutoff)

% inputs: 
    % r_0: vertical start postion (m) 3x1
    % r_t: target vertical end position (m) 3x1
    % v_0: vertical velocity start postion (m/s) 3x1
    % v_t: vertical velocity end position (m/s) 3x1
    % a_t: vertical target acceleration (m/s2) 3x1
    % t_go: time till target is reached (s)
    % t_cutoff: cutoff time to avoid dividing by vals close to zero (s)
    
 % output:
    % C_0: 3x1
    % C_1: 3x1
    % C_2: 3x1
    
    % to prevent dividing by near zero values of t_go as target is
    % approached
    if t_go < t_cutoff
        t_go = t_cutoff;
    end

    C_0 = a_t - (6 * ((v_t - v_0)/t_go) + (12 * (r_t - r_0) / t_go^2));
    C_1 = (-6 * a_t/t_go) + 6 * ((5 * v_t - 3 * v_0) / t_go^2) - (48 * ((r_t - r_0) / t_go^3));
    C_2 = (6 * a_t/t_go^2) - 12 * ((2 * v_t + v_0) / t_go^3) + (36 * ((r_t - r_0) / t_go^4));
    
end

function [a_prof, v_prof, r_prof] = getStateProfile(C_0, C_1, C_2, t, v_0, r_0)

% inputs: 
    % C_0: profile coefficient 3x1
    % C_1: profile coefficient 3x1
    % C_2: profile coefficient 3x1
    % t:   time since intial state (s)
    % r_0: vertical start postion (m) 3x1
    % v_0: vertical velocity start postion (m/s) 3x1
    
 % output:
    % a_prof: acceleration profile (m/s2) 3x1
    % v_prof: velocity profile (m/s)3x1
    % r_prof: position profile (m) 3x1
    
    a_prof = C_0 + (C_1 * t) + (C_2 * t^2); 
    v_prof = (C_0 * t) + (1/2 * C_1 * t^2) + (1/3 * C_2 * t^3) + v_0;
    r_prof = (1/2* C_0 * t^2) + (1/6 * C_1 * t^3) + (1/12 * C_2 * t^4) + (v_0 * t) + r_0;

end

