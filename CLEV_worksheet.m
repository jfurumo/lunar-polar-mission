%John Furumo
%ENAE 791
%Term Project
%Mass Estimation Worksheet

%clear;
clf;
clc;
close all;
format shortg;
set(0,'DefaultFigureWindowStyle','docked');

%% CLEV spacecraft parameters
mass=15000; %kg
diameter=4; %m
cD=1.2;
L_D=0.25; %lift to drag ratio (dimensionless)
%ballistic coefficient
beta=mass/(cD*pi*diameter^2/4); %kg/m^2
disp(strcat('ballistic coefficient, beta =',{' '},num2str(beta),{' '},'kg/m^2'))


%Orion lifting entry
P=2*pi*sqrt(a^3/mu_earth); %s
period=P/16; %s
gamma=deg2rad(FPA); %rad
state_ei=[r_ei;v_ei;gamma;theta;0;0;0];
L_D=0.25; %lift to drag ratio (dimensionless)
tspan=[0:1:5800]; %s
tolerance=1e-013;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
%2BP orbit propagation (numerical integration of 2BP), canonical coordinates
[t_lifting,state_lifting]=ode45(@integrate_2BP_canonical_LD,tspan,state_ei,options,mu_earth,g0,r_earth,rho_0,h_s,beta,L_D,state_ei(4));
%convert spacecraft position to altitude above Earth surface
state_lifting(:,1)=state_lifting(:,1)-r_earth;
%numerically differentiate velocity to get deceleration
for j=1:length(t_lifting)-1
    state_lifting(j,8)=abs(1000*(state_lifting(j+1,2)-state_lifting(j,2))); %m/s^2
end
state_lifting(end,8)=state_lifting(end-1,8);
peak_decel=max(state_lifting(:,8));
disp(strcat('peak deceleration =',{' '},num2str(peak_decel),{' '},'m/s^2'))
%terminal velocity
vT=sqrt(-(2*g0*(beta)*sin(-pi/2))/(rho_0)); %m/s
disp(strcat('terminal velocity =',{' '},num2str(vT),{' '},'m/s'))

%part a) Plot altitude v velocity
disp('part a)')
figure('name','Orion EFT-1 Lifting Entry')
subplot(2,2,1)
plot(state_lifting(:,2),state_lifting(:,1),'LineWidth',2)
xlabel('Velocity (km/s)')
ylabel('Altitude (km)')
ylim([0,130])
title('5a)')

%part b) Plot altitude v time
disp('part b)')
subplot(2,2,2)
plot(t_lifting,state_lifting(:,1),'LineWidth',2)
xlabel('Time (s)')
ylabel('Altitude (km)')
ylim([0,130])
title('5b)')

%part c) Plot altitude v downrange distance
disp('part c)')
subplot(2,2,3)
plot(state_lifting(:,7),state_lifting(:,1),'LineWidth',2)
xlabel('Downrange Distance (km)')
ylabel('Altitude (km)')
ylim([0,130])
title('5c)')

%part d) Plot deceleration v altitude
disp('part d)')
subplot(2,2,4)
plot(state_lifting(:,1),state_lifting(:,8),'LineWidth',2)
xlabel('Altitude (km)')
ylabel('Deceleration (m/s^2)')
title('5d)')