%John Furumo
%ENAE 601 section 0101
%Assignment 4
%due 11/25/19

clear;
clf;
clc;
close all;
format shortg;

mu_earth=398600.4415; %km^3/s^2
radius_earth=6378; %km
mu=mu_earth;

%Question 1: Earth-Moon 3BP
disp('Question 1')
%JPL Horizons data
m_earth=5.97219e+24; %kg
m_moon=7.349e+22; %kg
mu_calc=m_moon/(m_earth+m_moon)
mu_accept=0.012277471;

%%
%Question 2: Earth-Moon 3BP with spacecraft
disp('Question 2')
%units
DU=384400; %km
SU=1.023; %km/s
TU=4.348*86400; %s
%spacecraft initial state

% x0=1.061692; %DU
% y0=0; %DU
% z0=0; %DU
% xdot0=0; %SU
% ydot0=0.403877; %SU
% zdot0=0; %SU

v_GEO=3; %SU
v_LEO=6.9; %SU
dv=0;

x0=0.026; %DU
y0=0; %LEO radius (DU)
z0=0; %DU
xdot0=0;
ydot0=6.17;
zdot0=0; %SU

state0=[x0;y0;z0;xdot0;ydot0;zdot0];
%Numerical integration of 3BP
t_end=1; %TU
t_step=0.0001; %TU
tspan=[0:t_step:t_end]; %TU
tolerance=1e-013;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
[t,rv]=ode45(@integrate_3BP_ND,tspan,state0,options,mu_accept);
%convert to physical units
t=t*TU; %s
rv(:,1:3)=rv(:,1:3)*DU; %km
rv(:,4:6)=rv(:,4:6)*SU; %km/s
%plot spacecraft trajectory in XY plane
figure('Name','3BP 3D')
plot3(rv(:,1),rv(:,2),rv(:,3),'r-','linewidth',2)
hold on;
axis equal;
plot3(0,0,0,'bo','linewidth',2);
plot3(1*DU,0,0,'ko','linewidth',2);
plot3(rv(1,1),rv(1,2),rv(1,3),'g+','linewidth',2)
plot3(rv(end,1),rv(end,2),rv(end,3),'r+','linewidth',2)
title({'3 Body Problem','Spacecraft Trajectory through Earth-Moon System (XY Plane)'})
legend([{'Spacecraft Trajectory'},{'Earth'},{'Moon'},{'Starting Point'},{'Ending Point'}])
xlabel('X: Distance from Earth (km)')
ylabel('Y: Distance from Earth (km)')
zlabel('Z: Distance from Earth (km)')
% xlim([-1.0 1.25])
% ylim([-1.0 1.0])
az = 0;
el = 90;
view(az, el);
%output = integrate_3BP(t,state0,mu_accept);


