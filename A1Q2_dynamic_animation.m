%%%%%%%%
%{
Name: Katie Alderton
Date:27/03/2024

Input(s):
- angle theta [degrees]
- given values (lengths, m_arm, I_o, r_COM, etc.)
- angular velocity (w) and angular acceleration (a)
- F_muscle, R_par, R_perp

Output(s):
- Animation of bicep curl, showing the moving arm segment OAB and the
contraction/extention of the attached bicep muscle
- Organised Subplot of animated lines representing:
    - theta [degrees] vs t(s)
    - R_par [N] vs t (s)
    - R_perp [N] vs t (s)
    - F_muscle [N] vs t (s) 
%}
clear
clc
close all
%
% Prompt user for inputs for the angle theta (range of motion, degrees),
% number of reps per minute (min^-1) and the mass of the carried weight
% (m_weight).
theta_o_deg=input('Please input an angle theta in degrees, characterising the range of motion.\n');
while theta_o_deg<0||theta_o_deg>70
    fprintf('ERROR: That value for theta is invalid. Please respect the demands of the question.\n');
    theta_o_deg=input('Please input a valid angle for theta that is between 0 to 70 degrees.\n');
end
% The number of reps per minute
rep_permin=input('Please input the number of bicep curl repetitions completed per minute [min^-1].\n');
while rep_permin<=0  %later on in T calculations, we divide by rep_permin so it cannot be 0. 
    fprintf('ERROR: You can not complete a negative (or 0) number of reps. Please respect the demands of the question.\n');
    rep_permin=input('Please input a valid number for the number of reps completed per minute.\n');
end
% mass of carried weight
m_weight=input('Please input the mass of the carried weight, [kg].\n');
while m_weight<=0
    fprintf('ERROR: Mass can neither be 0 (as you would hence not be carrying any weight) or negative in the nature of this question.\n');
    m_weight=input('Please input a valid mass of the carried weight, [kg].\n');
end
%
% Define variables 
a=0.15; % length from O to A in metres
b=0.45; % length from O to B in metres
mu_u=0.20; % upper arm muscle attachment distance [m]
mu_d=0.03; %forearm muscle attachment distance [m]
m_arm=3.5; % arm mass [kg]
g=9.81; % gravitational acceleration on earth [m/s^2]
I_o=0.08875+0.2065*m_weight; % mass moment of inertia of the system about point O
r_COM=(0.525+0.45*m_weight)/(3.5+m_weight); % location of the centre of mass of the system
%
%convert theta_O [degs] to radians
theta_o_rad=deg2rad(theta_o_deg);
%
%set up t values
t=linspace(0,30,1000); %sets t values for the domain 0<=t<=30 [s].
%
% Angular displacement, theta(t)
T=60/(rep_permin); %T [s] is the period of motion (time for one complete cycle) 
phi=2*pi/T; % [per sec] is the angular frequency of the motion
theta_t=theta_o_rad*cos(phi*t); % gives the angular displacement in radians
theta_t_deg=rad2deg(theta_t); % angular displacement in degrees
%
% Angular velocity (w) and angular acceleration (ac) (will be vectors) 
omega=-theta_o_rad*phi*sin(phi*t); % angular velocity in [rad/s]
omega_deg=rad2deg(omega); % angular velocity [deg/s]
%
alpha=-theta_o_rad*(phi^2)*cos(phi*t); % angular acceleration in [rad/s^2]
alpha_deg=rad2deg(omega); % angular acceleration [deg/s^2]
%
% Compute the forces below (F_muscle, R_par, R_perp) by splitting them up
% to simplify and reduce errors. 
%
% Muscle Force, F_muscle
F_mus_numerator=sqrt((mu_u^2)+(mu_d^2)-2*mu_u*mu_d*sin(theta_t));
F_mus_denominator=mu_u*mu_d*cos(theta_t);
F_mus_bracket=I_o*alpha+(a*m_arm+b*m_weight)*g*cos(theta_t);
F_muscle=(F_mus_numerator./F_mus_denominator).*F_mus_bracket; %computes muscle force [N]
%
% Parallel reaction force 
R_par_numerator=mu_d-mu_u*sin(theta_t);
R_par_denominator=F_mus_numerator;
R_par_bracket=(m_arm+m_weight)*(g*sin(theta_t)-(omega.^2*r_COM));
R_par=F_muscle.*(R_par_numerator./R_par_denominator)+R_par_bracket; % parallel reaction force [N]
%
% Perpendicular reaction force 
R_perp_bracket=(m_arm+m_weight)*(g*cos(theta_t)+alpha*r_COM);
R_perp_numerator=mu_u*cos(theta_t);
R_perp_denominator=F_mus_numerator;
R_perp=R_perp_bracket-F_muscle.*(R_perp_numerator./R_perp_denominator);
% 
% define position of weight
weight_x=b*cos(theta_t); % defines x position of the carried weight at time t
weight_y=b*sin(theta_t); % defines y position of the carried weight at time t
%
% Muscle line position
muscle_x=0.03*cos(theta_t);
muscle_y=0.03*sin(theta_t);
%
% Set up Subplotting scheme 
% Bicep curl animation
subplot(4,2,[1,3,5,7]) % for bicep curl animation (LHS of 4x2 plot)
arm_OAB=animatedline('Color','#FCEACC','LineWidth',10);
muscle_line=animatedline('Color','r','LineWidth',4);
weightball=animatedline('Marker','o','MarkerSize',30,'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1.5);
title('Bicep Curl Animation')
axis equal
axis([0 b -b b])
%
% Animated plot of theta [degrees] vs. t [s]
subplot(4,2,2)
theta_line=animatedline('Color','b','LineWidth',1.5);
title('Graph of θ [deg] vs. t [s]')
xlabel('t [s]')
ylabel('θ [deg]')
axis([0 30 -theta_o_deg-0.1 theta_o_deg+0.1]) % the +-0.1 in the axes limit is set for the case that theta_o_rad is inputted as 0
%
% Animated plot of the parallel reaction force [N] vs. time [s]
subplot(4,2,4)
Rpar_line=animatedline('Color','r','LineWidth',1.5);
title('Graph of R_{par} [N] vs. t [s]')
xlabel('t [s]')
ylabel('R_{par} [N]')
axis([0 30 min(R_par)-0.1 max(R_par)+0.1]) % the +-0.1 in the axes limit is set for the case that theta_o_rad is inputted as 0
%
% Animated plot of the perpendicular reaction force [N] vs. time [s]
subplot(4,2,6)
Rperp_line=animatedline('Color','g','LineWidth',1.5);
title('Graph of R_{perp} [N] vs. t [s]')
xlabel('t [s]')
ylabel('R_{perp} [N]')
axis([0 30 min(R_perp)-0.1 max(R_perp)+0.1]) % the +-0.1 in the axes limit is set for the case that theta_o_rad is inputted as 0
%
% Animated plot of Muscle force [N] vs. t [s]
subplot(4,2,8)
Fmuscle_line=animatedline('Color','y','LineWidth',1.5);
title('Graph of F_{muscle} [N] vs. t [s]')
xlabel('t [s]')
ylabel('F_{muscle} [N]')
axis([0 30 min(F_muscle)-0.1 max(F_muscle)+0.1]) % the +-0.1 in the axes limit is set for the case that theta_o_rad is inputted as 0
%
%
n=length(t); % number of points in t
for i=1:n
    subplot(4,2,[1,3,5,7])
    clearpoints(muscle_line)
    clearpoints(arm_OAB) % clears previous points on the arm segment OAB
    clearpoints(weightball) % clears previous points in carried weight
    addpoints(arm_OAB ,[0,weight_x(i)],[0,weight_y(i)]); % adds points (2 points that connect a line) to arm segment OAB
    addpoints(muscle_line,[0,muscle_x(i)], [0.20,muscle_y(i)]) %adds points with a fixed position at (0,0.20)
    addpoints(weightball,weight_x(i),weight_y(i)) % adds points (only x-y position of the weight is relevant) to carried weight.
    drawnow
    %
    % Animates plot of angular displacement
    subplot(4,2,2)
    addpoints(theta_line,t(i),theta_t_deg(i))
    drawnow
    %
    % Animates plot of R_par [N] vs t [s]
    subplot(4,2,4)
    addpoints(Rpar_line,t(i),R_par(i))
    drawnow
    %
    % Animates plot of R_perp [N] vs t [s]
    subplot(4,2,6)
    addpoints(Rperp_line,t(i),R_perp(i))
    drawnow
    %
    % Animates plot of Muscle force [N] vs t [s]
    subplot(4,2,8)
    addpoints(Fmuscle_line,t(i),F_muscle(i))
    drawnow
end


