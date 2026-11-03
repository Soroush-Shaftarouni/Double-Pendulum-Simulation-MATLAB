clc; clear; close all;

p.m1 = 1; p.m2 = 1; p.L1 = 1; p.L2 = 1; p.g = 9.81;
z0 = [pi/2; 0; pi/2; 0]; 
tbaze = linspace(0, 20, 2001); 

[t, z] = ode45(@(t, z) pendulum_equations(t, z, p), tbaze, z0);

teta1 = z(:, 1); omega1 = z(:, 2);
teta2 = z(:, 3); omega2 = z(:, 4);

T = 0.5*(p.m1+p.m2)*p.L1^2.*omega1.^2 + 0.5*p.m2*p.L2^2.*omega2.^2 + p.m2*p.L1*p.L2.*omega1.*omega2.*cos(teta1-teta2);
V = -(p.m1+p.m2)*p.g*p.L1.*cos(teta1) - p.m2*p.g*p.L2.*cos(teta2);
E = T + V;

dt = tbaze(2) - tbaze(1);
sum_odd = sum(E(2:2:end-1));
sum_even = sum(E(3:2:end-2));
Integral_Simpson = (dt/3) * (E(1) + 4*sum_odd + 2*sum_even + E(end));
Average_E_Simpson = Integral_Simpson / tbaze(end);

Average_E_Matlab = mean(E);

fprintf('میانگین انرژی (سیمپسون): %f\n', Average_E_Simpson);
fprintf('میانگین انرژی (mean): %f\n', Average_E_Matlab);
