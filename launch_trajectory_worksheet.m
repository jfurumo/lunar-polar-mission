%John Furumo
%ENAE 791 section 0101
%Lunar Mission Project
%date 04/02/2020

clear;
clf;
clc;
close all;
%format shortg;
set(0,'DefaultFigureWindowStyle','docked');

%% CONSTANTS
G=6.67430e-20; %km^3/kg*s^2

%Earth
m_earth=5.97219e+24; %kg
mu_earth=398600.4415; %km^3/s^2
r_earth=6378; %km
v_equator=2*pi*r_earth/86400 %tangential Earth velocity at equator (km/s)
latitude=28.5; %deg
v_site=2*pi*(r_earth*cosd(latitude))/86400 %launch site tangential velocity due to Earth rotation (km/s)
g0_earth=G*m_earth/r_earth^2 %km/s^2

%Moon
m_moon=7.349e+22; %kg
mu_moon=4902.800066; %km^3/s^2 
r_moon=1737.4; %km
g0_moon=G*m_moon/r_moon^2 %km/s^2

%% Earth launch trajectory
g0=9.81; %m/s^2
r0=r_earth*1000; %m
rho_0=1.226; %kg/m^3
h_s=7524; %m
d_LV=4; %launch vehicle diameter (m)
A_LV=pi*(d_LV/2)^2; %launch vehicle cross-sectional area (m^2)
c_D=2; %drag coefficient
L_D=2; %lift-to-drag ratio
T=110100; %rocket engine thrust (N)
v_e=4565; %m/s
m_dot_prop=3.5+20.6; %propellant flow rate (kg/s)

%vertical launch segment
v_0=1e-04; %initial velocity (m/s)
m_0=10000; %kg
X_0_vert=[r0;v_0;0;0;m_0] %LV initial state vector
t_0_vert=0; %s
t_f_vert=60; %s
t_span_vert=[t_0_vert:1:t_f_vert]; %s
tolerance=1e-006;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
launch_state = integrate_vertical_launch_state(t_span_vert,X_0_vert,mu_earth,g0,r0,rho_0,h_s,A_LV,c_D,T,m_dot_prop)
[t_vert,rv_vert]=ode45(@integrate_vertical_launch_state,t_span_vert,X_0_vert,options,mu_earth,g0,r0,rho_0,h_s,A_LV,c_D,T,m_dot_prop);

%planar launch segment
t_0_planar=t_f_vert; %s
t_f_planar=t_0_planar+360; %s
t_span_planar=[t_0_planar:1:t_f_planar]; %s
gamma_0=pi/2; %rad
theta_0=0; %rad
X_0_planar=rv_vert(end,1:5)'
X_0_planar(3)=gamma_0
phi_0=deg2rad(0); %rad
U_0=[T;phi_0;v_e];
launch_state = integrate_planar_launch_state(t_span_planar,X_0_planar,U_0,mu_earth,r0,rho_0,h_s,A_LV,c_D,L_D)
[t_planar,rv_planar]=ode45(@integrate_planar_launch_state,t_span_planar,X_0_planar,options,U_0,mu_earth,r0,rho_0,h_s,A_LV,c_D,L_D);

figure('name','Launch Altitude')
plot(t_vert,rv_vert(:,1)-r_earth,'bo')
hold on;
plot(t_planar,rv_planar(:,1)-r_earth,'bo')
xlabel('Time (s)')
ylabel('Altitude (m)')
hold off;

figure('name','Launch Velocity')
plot(t_vert,rv_vert(:,2),'ro')
hold on;
plot(t_planar,rv_planar(:,2),'ro')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
hold off;

figure('name','Launch Mass')
plot(t_vert,rv_vert(:,5),'gs')
hold on;
plot(t_planar,rv_planar(:,5),'gs')
xlabel('Time (s)')
ylabel('Launch Vehicle Mass (kg)')
hold off;

%% Lunar ascent trajectory
%launch_state = integrate_planar_launch_state(t,X,U,mu,r0,rho_0,h_s,A,c_D,L_D);
%X - column vector with initial state (r,v,gamma,theta,m)
%U - column vector with control inputs (T,phi,v_e)
tolerance=1e-013;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
burnout=600; %s
t_0_asc=0; %s
t_f_asc=burnout; %s
t_span_lunar_asc=[0:1:burnout];
rho_0=0;
A=0;
c_D=0;
L_D=0;
h_s=1;
r0=0;
X_0_lunar=[r_moon;1e-04;pi/2;0;4555]
% X=[0;0;r_moon]; %km
% V=[0;0;0.1]; %km/s
% state_cart=[X;V];
% state_oe=cart2oe(state_cart,mu_earth);
%state_oe(1:5);
phi_0=deg2rad(5); %rad
U_0_lunar=[T;phi_0;v_e];
launch_state = integrate_planar_launch_state(t_span_lunar_asc,X_0_lunar,U_0_lunar,mu_moon,r0,rho_0,h_s,A,c_D,L_D)
[t_asc,rv_asc]=ode45(@integrate_planar_launch_state,t_span_lunar_asc,X_0_lunar,options,U_0_lunar,mu_moon,r0,rho_0,h_s,A,c_D,L_D);
% for j=1:length(rv_asc)
%     state_oe_canon(j,:)=[state_oe(1:5)',rad2deg(rv_asc(j,4))];
%     state_cart_canon(j,:) = (oe2cart(state_oe_canon(j,:)'))';
% end

figure('name','Lunar Ascent Altitude')
plot(t_asc,rv_asc(:,1)-r_earth,'bo')
xlabel('Time (s)')
ylabel('Altitude (m)')
hold off;

figure('name','Lunar Ascent Velocity')
plot(t_asc,rv_asc(:,2),'ro')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
hold off;

figure('name','Lunar Ascent Mass')
plot(t_asc,rv_asc(:,3),'gs')
xlabel('Time (s)')
ylabel('Launch Vehicle Mass (kg)')
hold off;