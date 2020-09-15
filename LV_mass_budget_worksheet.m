%John Furumo
%ENAE 791
%Term Project
%Launch Vehicle Mass Estimation Worksheet

clear;
clf;
clc;
close all;
format shortg;
set(0,'DefaultFigureWindowStyle','docked');

% %folder pathway
% folder=('C:\Users\jfuru\Desktop\UMD\Spring 2020\ENAE 791\Project\');
% %file name
% file=('space bus.tif');
% %full pathway
% path=[folder file];
% I=imread(path);
% figure('name','SpaceBus Diagram')
% imshow(I);
% text(10,100,'\leftarrow MNTF')



%% CONSTANTS
R_moon=1737.4; %km
mu_moon=4902.800066; %km^3/s^2 
g0_moon=mu_moon/(R_moon^2); %km/s^2
g0=9.81e-03; %km/s^2
rho_LO2=1140; %density of LO2 (kg/m^3)
rho_LH2=71; %density of LH2 (kg/m^3)
R=8.31446261815324; %ideal gas constant
M_He=4.002602e-03; %Helium molar mass (kg/mol)
delta_v_FOS=1.1; %factor of safety for delta-v budget

%% RL10B-2 engine (Hydrolox)
Isp_RL10=465.5; %specific impulse (s)
T_RL10=110100; %thrust (N) (=24800lbf)
v_exit_RL10=4.565; %exit velocity (km/s)
P_chamber_RL10=4412000; %combustion chamber pressure (Pa)
m_dot_prop=3.5+20.6; %propellant mass flow rate (kg/s)
m_RL10=301.2; %dry mass (kg)
MR_RL10=5.88; %mixture ratio

%% RS-25 engine (Hydrolox)
Isp_RS25_vac=452.3; %specific impulse, vacuum (s)
Isp_RS25_SL=366; %specific impulse, sea level (s)
T_RS25_vac=2279000; %thrust, vacuum (N) (=512300lbf)
T_RS25_SL=1860000; %thrust, sea level (N) (=418000lbf)
v_exit_RS25_vac=4.436; %exit velocity, vacuum (km/s)
v_exit_RS25_SL=3.59; %exit velocity, sea level (km/s)
P_chamber_RS25=20640000; %combustion chamber pressure (Pa) (=2994psi)
m_RS25=3527; %dry mass (kg)
MR_RS25=6.03; %mixture ratio

%% RS-68A engine (Hydrolox)
Isp_RS68_vac=412; %specific impulse, vacuum (s)
Isp_RS68_SL=362; %specific impulse, sea level (s)
T_RS68_vac=3559716; %thrust, vacuum (N) (=800000lbf)
T_RS68_SL=3137000; %thrust, sea level (N) (=705000lbf)
v_exit_RS68_vac=4.04; %exit velocity, vacuum (km/s)
v_exit_RS68_SL=3.551; %exit velocity, sea level (km/s)
P_chamber_RS68=10260000; %combustion chamber pressure (Pa) (=1488psi)
m_RS68=14870; %dry mass (kg)
MR_RS68=5.97; %mixture ratio

%% RD-180 engine (Kerolox)
Isp_RD180_vac=338; %specific impulse, vacuum (s)
Isp_RD180_SL=311; %specific impulse, sea level (s)
T_RD180_vac=4150000; %thrust, vacuum (N) (=930000lbf)
T_RD180_SL=3830000; %thrust, sea level (N) (=860000lbf)
v_exit_RD180_vac=3.31; %exit velocity, vacuum (km/s)
v_exit_RD180_SL=3.05; %exit velocity, sea level (km/s)
P_chamber_RD180=26700000; %combustion chamber pressure (Pa) (=3870psi)
m_RD180=5480; %dry mass (kg)
MR_RD180=2.72; %mixture ratio

%% LR87-5 engine (Hypergolic)
Isp_LR87_vac=297; %specific impulse, vacuum (s)
Isp_LR87_SL=259; %specific impulse, sea level (s)
T_LR87_vac=1096800; %thrust, vacuum (N) (=246279lbf)
T_LR87_SL=956500; %thrust, sea level (N) (=214775lbf)
v_exit_LR87_vac=2.914; %exit velocity, vacuum (km/s)
v_exit_LR87_SL=2.541; %exit velocity, sea level (km/s)
P_chamber_LR87=5400000; %combustion chamber pressure (Pa) (=782.7psi)
m_LR87=739; %dry mass (kg)
MR_LR87=1.93; %mixture ratio


%% SSTO Earth Launch Vehicle Design
m_payload=13500; %launch vehicle payload to LEO (kg)
%m_payload=25000;
dv_LEO=9.3; %launch to LEO delta-v estimate (km/s)

%Hydrolox - 4x RS-25
MR=MR_RS25; %RS-25 mixture ratio
%m_LV_wet=m_LV_dry*exp((dv_launch)/(Isp*g0));
Isp_launch=(Isp_RS25_vac+Isp_RS25_SL)/2; %average of sea level and vacuum Isp values (s)
T_launch=(T_RS25_vac+T_RS25_SL)/2; %average of sea level and vacuum thrust values (N)
dv_launch=delta_v_FOS*dv_LEO; %delta-v for launch from Earth surface to LEO (km/s)
d_LV=4; %vehicle diameter (m)
n_engine=6; %number of engines
m_LV_engine=n_engine*m_RS25*1.5; %four RS25 engines plus associated plumbing and other hardware (kg)
%Helium pressurization system
n_LV_COPV=4; %Helium pressurization system COPV quantity, based on Centaur
d_COPV=0.66; %Helium pressurization system COPV diameter, based on Centaur (m)
V_LV_COPV=n_LV_COPV*(4/3)*pi*(d_COPV/2)^3; %Helium pressurization system COPV volume, based on Centaur (m^3)
m_LV_COPV=115.3*V_LV_COPV+3; %Helium pressurization system COPV mass, based on Centaur (kg)
p_COPV=4000*6894.76; %COPV Helium pressure (Pa)
T_COPV=275; %COPV Helium temperature (K)
m_LV_He=(p_COPV*V_LV_COPV*M_He)/(R*T_COPV); %mass of Helium pressurization gas (kg)
m_LV_COPV=m_LV_COPV+m_LV_He; %total pressurization system mass (kg)
%propellant tanks
m_LV_tank_estimate=2000; %first estimate for propellant tank mass (kg)
m_LV_tank=0; %calculated propellant tank mass (kg)
m_LV_dry=m_payload+m_LV_engine+m_LV_COPV+m_LV_tank_estimate; %dry mass of LV (kg)
tol=1e-03; %convergence tolerance
dm_tank=m_LV_tank-m_LV_tank_estimate; %tank mass error
count=0;
while abs(dm_tank)>tol
    %propellant estimation
        m_LV_wet=m_LV_dry*exp((dv_launch)/(Isp_launch*g0)); %LV total mass (kg)
        m_LV_prop=m_LV_wet-m_LV_dry; %LV propellant mass (kg)
    %tankage
        %oxidizer
        m_LV_LO2=MR*(m_LV_prop/(MR+1)); %mass of LO2 required (kg)
        V_LV_LO2=m_LV_LO2/rho_LO2; %LO2 tank volume (m^3)
        [A_tank_LV_LO2,l_tank_LV_LO2]=tank_area(V_LV_LO2,d_LV); %LO2 tank surface area (m^2)
        m_insulation_LV_LO2=1.123*A_tank_LV_LO2; %LO2 tank insulation mass (kg)    
        m_LV_tank_LO2=12.16*V_LV_LO2; %LO2 tank mass estimate (kg)
        %fuel
        m_LV_LH2=(m_LV_prop/(MR+1)); %mass of LH2 required (kg)
        V_LV_LH2=m_LV_LH2/rho_LH2; %LH2 tank volume (m^3)
        [A_tank_LV_LH2,l_tank_LV_LH2]=tank_area(V_LV_LH2,d_LV); %LH2 tank surface area (m^2)
        m_insulation_LV_LH2=2.88*A_tank_LV_LH2; %LH2 tank insulation mass (kg)
        m_LV_tank_LH2=9.09*V_LV_LH2; %LH2 tank mass estimate (kg)
        %total
        m_LV_tank=m_LV_tank_LO2+m_insulation_LV_LO2+m_LV_tank_LH2+m_insulation_LV_LH2; %total LV propellant tank mass (kg)
    %fairings
        A_fore=pi*d_LV*(d_LV/2); %foreward fairing area (m^2)
        A_intertank=pi*d_LV*(d_LV); %intertank fairing area (m^2)
        A_aft=pi*d_LV*(d_LV/2); %aft fairing area (m^2)
        A_LV_fairing=A_fore+A_intertank+A_aft; %total LV fairing surface area (m^2)
        m_LV_fairing=4.95*(A_LV_fairing)^1.15; %mass of LV fairing (kg)
    %avionics
        m_LV_avionics=10*m_LV_wet^0.361; %mass of avionics (kg)
    %wiring
        l_LV=l_tank_LV_LO2+l_tank_LV_LH2; %LV length (m)
        m_LV_wiring=1.058*sqrt(m_LV_wet)*l_LV^0.25; %mass of wiring (kg)
    %thrust structure
        m_LV_thrust=n_engine*(2.55e-04)*T_launch; %thrust structure mass (kg)
    
    dm_tank=m_LV_tank-m_LV_tank_estimate; %propellant tank mass change (kg)
    m_LV_tank_estimate=m_LV_tank; %new calculated mass becomes new mass estimate (kg)
    m_LV_dry=m_payload+m_LV_engine+m_LV_COPV+m_LV_tank+m_LV_fairing+m_LV_avionics+m_LV_wiring+m_LV_thrust;
    count=count+1;
end
disp(strcat('LV - SSTO Hydrolox,', num2str(n_engine), 'x RS-25'))
m_LV_dry %dry mass of Space Bus (kg)
m_LV_wet=m_LV_dry*exp((dv_launch)/(Isp_launch*g0)) %Space Bus wet mass (kg)
m_LV_prop=m_LV_wet-m_LV_dry %Space Bus propellant mass (kg)
TW_ratio_LV=(n_engine*T_RS25_SL)/(m_LV_wet*g0*1000) %LV thrust to weight ratio

%Kerolox

%Hypergolic 





%% FUNCTIONS

%this function calculates the surface area of a propellant tank based on its
%volume
%INPUTS
%V - volume (m^3)
%d - tank diameter (m)
function [A_tank,l_tank] = tank_area(V_total,d)
A_dome=4*pi*(d/2)^2; %dome surface area (m^2)
V_dome=(4/3)*pi*(d/2)^3; %dome volume (m^3)
V_cylinder=V_total-V_dome; %cylinder volume (m^3)
l_cylinder=V_cylinder/(pi*(d/2)^2); %cylinder length (m)
A_cylinder=pi*d*l_cylinder; %cylinder surface area (m^2)
A_tank=A_cylinder+A_dome; %total tank surface area (m^2)
l_tank=l_cylinder+d; %total tank length (m)
end

%% ARCHIVE

% Space Bus - Destiny ISS module
% dv_SB=9.56; %delta-v for TLI+TEI round trip (km/s)
% m_module=14.5*1000; %based on ISS Destiny module (kg)
% m_engine=m_RL10*1.5; %one RL10 engine plus associated plumbing and other hardware (kg)
% m_SB_tank=3673 %first estimate for propellant tank mass (kg)
% m_SB_dry=m_module+m_engine+m_SB_tank %dry mass of Space Bus (kg)
% m_SB_tot=m_SB_dry*exp((dv_SB)/(Isp*g0)); %Space Bus total mass (kg)
% m_SB_prop=m_SB_tot-m_SB_dry %Space Bus propellant mass (kg)
% 
% m_SB_LO2=MR*(m_SB_prop/(MR+1)); %mass of LO2 required (kg)
% m_SB_LH2=(m_SB_prop/(MR+1)); %mass of LH2 required (kg)
% 
% rho_LO2=1140; %density of LO2 (kg/m^3)
% rho_LH2=71; %density of LH2 (kg/m^3)
% 
% V_LO2=m_SB_LO2/rho_LO2; %LO2 tank volume (m^3)
% V_LH2=m_SB_LH2/rho_LH2; %LH2 tank volume (m^3)
% 
% m_SB_tank_LO2=12.16*V_LO2; %LO2 tank mass estimate (kg)
% m_SB_tank_LH2=9.09*V_LH2; %LH2 tank mass estimate (kg)
% 
% m_SB_tank_estimate=m_SB_tank_LO2+m_SB_tank_LH2 %total Space Bus propellant tank mass (kg)


% Delta-v budget
% %lunar ascent delta-v
% v_lunar_park=sqrt(mu_moon/r_park); %km/s
% dv_lunar_asc=v_lunar_park;
% %propellant mass required (work backwards through maneuvers)
% m_asc_tot=m_asc_dry*exp((dv_lunar_asc)/(Isp*g0)); %kg
% m_prop_asc=m_asc_tot-m_asc_dry; %kg
% 
% %lunar descent delta-v
% dv_lunar_desc=v_lunar_park;
% %propellant mass required (work backwards through maneuvers)
% m_desc_tot=(m_desc_dry+m_asc_dry+m_prop_asc)*exp((dv_lunar_desc)/(Isp*g0)); %kg
% m_prop_desc=m_desc_tot-(m_desc_dry+m_asc_dry+m_prop_asc); %kg
% 
% dv_TLI=5.5561; %km/s
% m_LM_tot=m_desc_tot*exp((dv_TLI)/(Isp*g0)); %kg
% m_prop_TLI=m_LM_tot-m_desc_tot; %kg