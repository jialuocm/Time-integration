% =========================================================================
% This is FEM Matlab code for dynamic analysis by direction intergration:
% (1) Generalized-alpha method, numerical dissipation controlled by user;
% (2) Two degree of freedom model ;  
% (3) Reference to J.Chung G.M. Hulbert paper 1993 ;  
% -------------------------------------------------------------------------
% By Jia Luo, 2023 May. 5th.
% =========================================================================
clear all; clc;

% Setup the alpha parameters;
% Central difference: beta = 0, gamma = 0.5
% Trapezoidal: beta = 0.25, gamma = 0.5
% Damped Newmark: beta = 0.3025, gamma = 0.6

rho_inf = 0.8;

alpha_m = (2 * rho_inf - 1) / (rho_inf + 1);
alpha_f = rho_inf / (rho_inf + 1);

gamma = 0.5 - alpha_m + alpha_f;
beta  = 0.25 *(1 - alpha_m + alpha_f)^2;


%beta  = 0.25 *(1 - alpha_m + alpha_f)^2;

% beta = 0;
% gamma = 0.5;

% beta = 0.25;
% gamma = 0.5;

% beta = 0.3025;
% gamma = 0.6;

m1 = 1;
m2 = 1;

k1 = 1.0e4;
k2 = 1.0;

% stiffness and mass matrix
M = [m1, 0.0; 0.0, m2];
K = [(k1 + k2), -k2; -k2, k2];

% initial condition
d0 = [ 1 ; 10];
v0 = [ 0 ; 0];

% determine the natrual frequencies
lambda = eig(K, M);
omega = sqrt(lambda);

T1 = 2 * pi / omega(1);
T2 = 2 * pi / omega(2);

%dt = T1 / 20;
%dt = T2 / pi ;
dt = 0.018;
T_period = 5 * T1;

N = ceil(T_period / dt);

%initial acceleration
a0 = M \ (-K * d0);

% allocate solutions
a = zeros(2, N+1);
v = a;
d = a;

a(:,1) = a0;
v(:,1) = v0;
d(:,1) = d0;
% LEFT, effecive mass matrix; 
LEFT = M * (1 - alpha_m) + (1 - alpha_f ) * beta * dt * dt * K;

for n = 2 : N + 1
    % Predictors, (9.1.7),(9.1.8) of linear FEM book;
    d_tilde = d(:,n-1) + dt * v(:,n-1) + 0.5 * dt * dt * (1 - 2 * beta) * a(:,n-1);
    v_tilde = v(:,n-1) + (1 - gamma) * dt * a(:,n-1);

    % RIGHT, effecive load at time tn;
    RIGHT = - M * alpha_m * a(:,n-1) - (1 - alpha_f ) * K * d_tilde - K * alpha_f * d(:,n-1);
    % LDL decomposition;
    a(:,n) = LEFT \ RIGHT;

    % Correctors;

    d(:,n) = d_tilde + beta * dt * dt * a(:,n);
    v(:,n) = v_tilde + gamma * dt * a(:,n);
end

% Postprocessing: Visualization

t = 1 : 1 : N;

subplot(2,2,1), plot(t, d(1,1:N)); grid on;
subplot(2,2,2), plot(t, v(1,1:N)); grid on;

subplot(2,2,3), plot(t, d(2,1:N)); grid on;
subplot(2,2,4), plot(t, v(2,1:N)); grid on;


%END