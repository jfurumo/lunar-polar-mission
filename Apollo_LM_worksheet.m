clear;
clf;
clc;
close all;
format shortg;
set(0,'DefaultFigureWindowStyle','docked');

%% CONSTANTS
r_moon=1737.4; %km
mu_moon=4902.800066; %km^3/s^2
g0=9.81e-03; %km/s^2

%% Apollo 11 Lunar Module
%descent stage
m_desc_dry=4483/2.20462; %kg
m_desc_prop=18184/2.20462; %kg
m_desc_tot=(m_desc_dry+m_desc_prop); %kg
Isp_desc=311; %specific impulse (s)
T_desc=45040; %thrust (N) = 10125 lbf
delta_v_desc_ref=2.5 %quoted delta-v for descent stage (km/s)
%ascent stage
m_asc_dry=4804/2.20462; %kg
m_asc_prop=5238/2.20462; %kg
m_asc_tot=(m_asc_dry+m_asc_prop); %kg
m_equip=569/2.20462; %kg
Isp_asc=311; %specific impulse (s)
T_asc=16000; %thrust (N) = 3500 lbf
delta_v_asc_ref=2.22 %quoted delta-v for descent stage (km/s)

m_total=m_desc_tot+m_asc_tot+m_equip; %LM total mass (kg)

%% Delta-v budget
%lunar ascent delta-v
delta_v_asc=-Isp_asc*g0*log(m_asc_dry/m_asc_tot)
%lunar descent delta-v
delta_v_desc=-Isp_desc*g0*log((m_asc_tot+m_desc_dry)/m_total)

% %propellant mass required (work backwards through maneuvers)
% m_desc_tot=(m_desc_dry+m_asc_dry+m_prop_asc)*exp((dv_lunar_desc)/(Isp*g0)); %kg
% m_prop_desc=m_desc_tot-(m_desc_dry+m_asc_dry+m_prop_asc); %kg
% 
% dv_TLI=5.5561; %km/s
% m_LM_tot=m_desc_tot*exp((dv_TLI)/(Isp*g0)); %kg
% m_prop_TLI=m_LM_tot-m_desc_tot; %kg