%%%%%%%%
%{
Name: Katie Alderton
Date:24/03/2024

Input(s):
- angle theta 
- given values for m values and lengths a and b

Output(s):
- Parallel Reaction force at O, (2 d.p)
- Perpendicular reaction force at O (2 d.p)
- Moment at O (2 d.p)
%}
clc
clear
close all
%
% Ask user to input an angle theta ("arbitrary" with respect to horizontal
% axis, so no angular limit is neeeded to be set by bounds of a typical bicep curl motion) 
theta_deg=input('Please enter an angle theta in [degrees].\n');

% convert this angle to radians
theta_rad=deg2rad(theta_deg);

% Define given values
m_arm=3.5; %[kg]
m_weight=20; %[kg]
a=15/100; %[m]
b=45/100; %[m]
g=9.81; %[m/s^2]

% Calculate Weight forces
W_arm=m_arm*g; %[N]
W_weight=m_weight*g; %[N]

% Calculate reaction forces and moment at point O
R_o_parallel=W_arm*sin(theta_rad)+W_weight*sin(theta_rad);
R_o_perpendicular=W_arm*cos(theta_rad)+W_weight*cos(theta_rad);
M_o=a*W_arm*cos(theta_rad)+b*W_weight*cos(theta_rad);

% Display results accurate to 2 decimal places
fprintf('Reaction force at point O parallel to the arm segment:\n %.2f [N].\n',R_o_parallel);
fprintf('Reaction force at point O perpendicular to the arm segment:\n %.2f [N].\n',R_o_perpendicular);
fprintf('Moment at point O:\n %.2f [Nm].\n',M_o);
