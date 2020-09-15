%John Furumo
%ENAE 791
%Term Project
%Mass Estimation Worksheet

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

%% RL10B-2 engine
Isp=465.5; %specific impulse (s)
T=110100; %thrust (N) (=24800lbf)
P_chamber=4412000; %combustion chamber pressure (Pa)
v_e=4.565; %exit velocity (km/s)
m_dot_prop=3.5+20.6; %propellant mass flow rate (kg/s)
m_RL10=301.2; %dry mass (kg)
MR=5.88; %mixture ratio

%% Lunar Scouting Rovers
m_Curiosity=1000; %lunar rover mass (kg)
m_Lunokhod=750; %lunar rover mass (kg)
m_Change4=1200; %lunar lander + rover mass (kg)

%% Space Bus - Destiny ISS Module
dv_bus=9.56; %delta-v for TLI+TEI round trip (km/s)
d_bus=4; %vehicle diameter (m)
V_destiny=104.77; %pressurized volume (m^3)
m_bus_crew=14.5*1000; %based on ISS Destiny module (kg)
m_engine=m_RL10*1.5; %one RL10 engine plus associated plumbing and other hardware (kg)
m_bus_tank_estimate=5000; %first estimate for propellant tank mass (kg)
m_bus_tank=0; %calculated propellant tank mass (kg)
tol=1e-03; %convergence tolerance
dm_tank=m_bus_tank-m_bus_tank_estimate; %tank mass error
count=0;
while abs(dm_tank)>tol
    m_bus_dry=m_bus_crew+m_engine+m_bus_tank_estimate; %dry mass of Space Bus (kg)
    m_bus_wet=m_bus_dry*exp((dv_bus)/(Isp*g0)); %Space Bus total mass (kg)
    m_bus_prop=m_bus_wet-m_bus_dry; %Space Bus propellant mass (kg)
    
    %oxidizer
    m_bus_LO2=MR*(m_bus_prop/(MR+1)); %mass of LO2 required (kg)
    V_LO2=m_bus_LO2/rho_LO2; %LO2 tank volume (m^3)
    [A_tank_bus_LO2,l_tank_bus_LO2]=tank_area(V_LO2,d_bus); %LO2 tank surface area (m^2)
    m_insulation_bus_LO2=1.123*A_tank_bus_LO2; %LO2 tank insulation mass (kg)    
    m_bus_tank_LO2=12.16*V_LO2; %LO2 tank mass estimate (kg)
    
    %fuel
    m_bus_LH2=(m_bus_prop/(MR+1)); %mass of LH2 required (kg)
    V_LH2=m_bus_LH2/rho_LH2; %LH2 tank volume (m^3)
    [A_tank_bus_LH2,l_tank_bus_LH2]=tank_area(V_LH2,d_bus); %LH2 tank surface area (m^2)
    m_insulation_bus_LH2=2.88*A_tank_bus_LH2; %LH2 tank insulation mass (kg)
    m_bus_tank_LH2=9.09*V_LH2; %LH2 tank mass estimate (kg)
    
    l_bus=l_tank_bus_LO2+l_tank_bus_LH2 %SpaceBus length (m)
    m_bus_tank=m_bus_tank_LO2+m_insulation_bus_LO2+m_bus_tank_LH2+m_insulation_bus_LH2; %total Space Bus propellant tank mass (kg)
    dm_tank=m_bus_tank-m_bus_tank_estimate; %propellant tank mass change (kg)
    m_bus_tank_estimate=m_bus_tank; %new calculated mass becomes new mass estimate (kg)
    count=count+1;
end
disp('Space Bus - Destiny')
m_bus_dry=m_bus_crew+m_engine+m_bus_tank %dry mass of Space Bus (kg)
m_bus_wet=m_bus_dry*exp((dv_bus)/(Isp*g0)) %Space Bus total mass (kg)
m_bus_prop=m_bus_wet-m_bus_dry %Space Bus propellant mass (kg)

%% Space Bus - Crew Transport (NASA JSC-26098)
dv_bus=delta_v_FOS*9.56; %delta-v for TLI+TEI round trip (km/s)
dv_TLI=delta_v_FOS*6.35; %delta-v from LEO to LLO
d_bus=4; %vehicle diameter (m)
n_engine=1; %number of engines
m_bus_engine=n_engine*m_RL10*1.5; %one RL10 engine plus associated plumbing and other hardware (kg)
n_bus_COPV=4; %Helium pressurization system COPV quantity, based on Centaur
d_COPV=0.66; %Helium pressurization system COPV diameter, based on Centaur (m)
V_bus_COPV=n_bus_COPV*(4/3)*pi*(d_COPV/2)^3; %Helium pressurization system COPV volume, based on Centaur (m^3)
m_bus_COPV=115.3*V_bus_COPV+3; %Helium pressurization system COPV mass, based on Centaur (kg)
p_COPV=4000*6894.76; %COPV Helium pressure (Pa)
T_COPV=275; %COPV Helium temperature (K)
m_bus_He=(p_COPV*V_bus_COPV*M_He)/(R*T_COPV); %mass of Helium pressurization gas (kg)
m_bus_COPV=m_bus_COPV+m_bus_He; %total pressurization system mass (kg)
%crew module
C_bus=4; %crew size (number of astronauts)
D_bus=7; %mission duration (days)
V_bus_crew=1.74*C_bus*D_bus^0.7444; %required volume for space bus crew module (m^3)
l_bus_crew=V_bus_crew/(pi*(d_bus/2)^2); %crew module length (m)
m_bus_crew=460*V_bus_crew^0.76; %mass of spacebus crew module (kg)

m_bus_tank_estimate=2000; %first estimate for propellant tank mass (kg)
m_bus_tank=0; %calculated propellant tank mass (kg)
m_bus_dry=m_bus_crew+m_bus_engine+m_bus_COPV+m_bus_tank_estimate; %dry mass of Space Bus (kg)
tol=1e-03; %convergence tolerance
dm_tank=m_bus_tank-m_bus_tank_estimate; %tank mass error
count=0;
while abs(dm_tank)>tol
    %propellant estimation
        m_bus_wet=m_bus_dry*exp((dv_bus)/(Isp*g0)); %Space Bus total mass (kg)
        m_bus_prop=m_bus_wet-m_bus_dry; %Space Bus propellant mass (kg)
    %tankage
        %oxidizer
        m_bus_LO2=MR*(m_bus_prop/(MR+1)); %mass of LO2 required (kg)
        V_bus_LO2=m_bus_LO2/rho_LO2; %LO2 tank volume (m^3)
        [A_tank_bus_LO2,l_tank_bus_LO2]=tank_area(V_bus_LO2,d_bus); %LO2 tank surface area (m^2)
        m_insulation_bus_LO2=1.123*A_tank_bus_LO2; %LO2 tank insulation mass (kg)    
        m_bus_tank_LO2=12.16*V_bus_LO2; %LO2 tank mass estimate (kg)
        %fuel
        m_bus_LH2=(m_bus_prop/(MR+1)); %mass of LH2 required (kg)
        V_bus_LH2=m_bus_LH2/rho_LH2; %LH2 tank volume (m^3)
        [A_tank_bus_LH2,l_tank_bus_LH2]=tank_area(V_bus_LH2,d_bus); %LH2 tank surface area (m^2)
        m_insulation_bus_LH2=2.88*A_tank_bus_LH2; %LH2 tank insulation mass (kg)
        m_bus_tank_LH2=9.09*V_bus_LH2; %LH2 tank mass estimate (kg)
        %total
        m_bus_tank=m_bus_tank_LO2+m_insulation_bus_LO2+m_bus_tank_LH2+m_insulation_bus_LH2; %total Space Bus propellant tank mass (kg)
    %fairings
        A_fore=pi*d_bus*(d_bus/2); %foreward fairing area (m^2)
        A_intertank=pi*d_bus*(d_bus); %intertank fairing area (m^2)
        A_aft=pi*d_bus*(d_bus/2); %aft fairing area (m^2)
        A_bus_fairing=A_fore+A_intertank+A_aft; %total bus fairing surface area (m^2)
        m_bus_fairing=4.95*(A_bus_fairing)^1.15; %mass of bus fairing (kg)
    %avionics
        m_bus_avionics=10*m_bus_wet^0.361; %mass of avionics (kg)
    %wiring
        l_bus=l_tank_bus_LO2+l_tank_bus_LH2; %SpaceBus length (m)
        m_bus_wiring=1.058*sqrt(m_bus_wet)*l_bus^0.25; %mass of wiring (kg)
    %thrust structure
        m_bus_thrust=n_engine*(2.55e-04)*T; %thrust structure mass (kg)
    
    dm_tank=m_bus_tank-m_bus_tank_estimate; %propellant tank mass change (kg)
    m_bus_tank_estimate=m_bus_tank; %new calculated mass becomes new mass estimate (kg)
    m_bus_dry=m_bus_crew+m_bus_engine+m_bus_COPV+m_bus_tank+m_bus_fairing+m_bus_avionics+m_bus_wiring+m_bus_thrust;
    count=count+1;
end
disp('Space Bus - JSC-26098')
m_bus_dry %dry mass of Space Bus (kg)
m_bus_wet=m_bus_dry*exp((dv_bus)/(Isp*g0)) %Space Bus wet mass (kg)
m_bus_prop=m_bus_wet-m_bus_dry %Space Bus propellant mass (kg)

%% Space Bus - Tanker (NASA JSC-26098)
m_tanker_dry=m_bus_dry-m_bus_crew-m_bus_engine-m_bus_thrust; %tanker dry mass (kg)
%m_tanker_transport=m_bus_launch-m_bus_tanker_dry %tanker propellant transfer mass (kg)
m_tanker_wet=20000; %launch vehicle payload mass to LEO (kg)
m_tanker_prop=m_tanker_wet-m_tanker_dry; %residual propellant to be carried in tanker (kg)
m_tanker_LO2=MR*(m_tanker_prop/(MR+1)); %mass of LO2 (kg)
m_tanker_LH2=(m_tanker_prop/(MR+1)); %mass of LH2 (kg)

m_bus_prop_residual=m_tanker_wet-m_bus_dry; %residual propellant mass in SpaceBus at launch (kg)

%Space Bus initial launch
disp('SpaceBus - One-way')
m_bus_launch=(m_bus_dry)*exp((dv_TLI)/(Isp*g0)) %SpaceBus launch mass (kg)
m_bus_launch_prop=m_bus_launch-m_bus_dry %SpaceBus launch mass, propellant (kg)
n_tanker_launches=ceil((m_bus_launch_prop-m_bus_prop_residual)/m_tanker_prop) %number of tanker flights needed

disp('SpaceBus - Round-trip')
m_bus_launch=(m_bus_dry)*exp((dv_bus)/(Isp*g0)) %SpaceBus launch mass (kg)
m_bus_launch_prop=m_bus_launch-m_bus_dry %SpaceBus launch mass, propellant (kg)
n_tanker_launches=ceil((m_bus_launch_prop-m_bus_prop_residual)/m_tanker_prop) %number of tanker flights needed

%% SSTO Lunar Descent/Ascent Vehicle
dv_LLV_desc=delta_v_FOS*1.6335; %delta-v to lunar surface from polar LLO (km/s)
dv_LLV_asc=delta_v_FOS*1.6335; %delta-v to reach polar LLO from lunar surface (km/s)
d_LLV=4; %vehicle diameter (m)
n_engine=4; %number of engines
m_LLV_engine=n_engine*m_RL10*1.5; %two RL10 engines plus associated plumbing and other hardware (kg)
n_LLV_COPV=4; %Helium pressurization system COPV quantity, based on Centaur
d_COPV=0.66; %Helium pressurization system COPV diameter, based on Centaur (m)
V_LLV_COPV=n_LLV_COPV*(4/3)*pi*(d_COPV/2)^3; %Helium pressurization system COPV volume, based on Centaur (m^3)
m_LLV_COPV=115.3*V_LLV_COPV+3; %Helium pressurization system COPV mass, based on Centaur (kg)
p_COPV=4000*6894.76; %COPV Helium pressure (Pa)
T_COPV=275; %COPV Helium temperature (K)
m_LLV_He=(p_COPV*V_LLV_COPV*M_He)/(R*T_COPV); %mass of Helium pressurization gas (kg)
m_LLV_COPV=m_LLV_COPV+m_LLV_He; %total pressurization system mass (kg)
%crew module
C_LLV=4; %crew size (number of astronauts)
D_LLV=14; %mission duration (days)
V_LLV_crew=1.74*C_LLV*D_LLV^0.7444; %required volume for LLV crew module (m^3)
l_LLV_crew=V_LLV_crew/(pi*(d_LLV/2)^2); %crew module length (m)
m_LLV_crew=460*V_LLV_crew^0.76; %mass of LLV crew module (kg) 

m_LLV_tank_LO2_estimate=5000; %first estimate for LO2 tank mass (kg)
m_LLV_tank_LO2=0; %calculated LO2 tank mass (kg)

m_LLV_tank_LH2_estimate=5000; %first estimate for LH2 tank mass (kg)
m_LLV_tank_LH2=0; %calculated LH2 tank mass (kg)
m_legs=1000; %mass of landing gear estimate (kg)
m_LLV_dry=m_LLV_crew+m_LLV_engine+m_LLV_COPV+m_legs+m_bus_tank+m_LLV_tank_LO2_estimate+m_LLV_tank_LH2_estimate; %LLV dry mass (kg)
% m_LLV_tank_estimate=m_LLV_tank_desc_estimate+m_LLV_tank_asc_estimate;
% m_LLV_tank=m_LLV_tank_desc+m_LLV_tank_asc;
tol=1e-03; %convergence tolerance
dm_tank_LO2=m_LLV_tank_LO2-m_LLV_tank_LO2_estimate;
dm_tank_LH2=m_LLV_tank_LH2-m_LLV_tank_LH2_estimate;
count=0;
while abs(dm_tank_LO2)>tol & abs(dm_tank_LH2)>tol
    %propellant estimation  
        %descent mass calculations
        m_LLV_desc=m_LLV_dry*exp((dv_LLV_desc)/(Isp*g0)); %LLV total mass at prior to descent (kg)
        m_LLV_prop_desc=m_LLV_desc-m_LLV_dry; %descent propellant mass (kg)
        %oxidizer
        m_LLV_LO2_desc=MR*(m_LLV_prop_desc/(MR+1)); %mass of LO2 required (kg)
        V_LLV_LO2_desc=m_LLV_LO2_desc/rho_LO2; %LO2 tank volume (m^3)
        m_LLV_tank_LO2_desc=12.16*V_LLV_LO2_desc; %LO2 tank mass estimate (kg)   
        %fuel
        m_LLV_LH2_desc=(m_LLV_prop_desc/(MR+1)); %mass of LH2 required (kg)
        V_LLV_LH2_desc=m_LLV_LH2_desc/rho_LH2; %LH2 tank volume (m^3)
        m_LLV_tank_LH2_desc=9.09*V_LLV_LH2_desc; %LH2 tank mass estimate (kg)
        %ascent mass calculations
        m_LLV_asc=(m_LLV_dry+m_bus_prop+m_LLV_prop_desc)*exp((dv_LLV_asc)/(Isp*g0)); %LLV total mass at liftoff (kg)
        m_LLV_prop_asc=m_LLV_asc-(m_LLV_dry+m_bus_prop+m_LLV_prop_desc); %ascent propellant mass (kg)
        %oxidizer
        m_LLV_LO2_asc=MR*(m_LLV_prop_asc/(MR+1)); %mass of LO2 required (kg)
        V_LLV_LO2_asc=m_LLV_LO2_asc/rho_LO2; %LO2 tank volume (m^3)
        m_LLV_tank_LO2_asc=12.16*V_LLV_LO2_asc; %LO2 tank mass estimate (kg)
        %fuel
        m_LLV_LH2_asc=(m_LLV_prop_asc/(MR+1)); %mass of LH2 required (kg)
        V_LLV_LH2_asc=m_LLV_LH2_asc/rho_LH2; %LH2 tank volume (m^3)
        m_LLV_tank_LH2_asc=9.09*V_LLV_LH2_asc; %LH2 tank mass estimate (kg)
        %propellant tank insulation
        %oxidizer
        V_LLV_LO2=V_LLV_LO2_desc+V_LLV_LO2_asc+V_bus_LO2; %total LO2 volume for ascent/descent + transfer (m^3)
        [A_tank_LLV_LO2,l_tank_LLV_LO2]=tank_area(V_LLV_LO2,d_LLV); %LO2 tank surface area (m^2)
        m_insulation_LLV_LO2=1.123*A_tank_LLV_LO2; %LO2 tank insulation mass (kg)    
        %fuel
        V_LLV_LH2=V_LLV_LH2_desc+V_LLV_LH2_asc+V_bus_LH2; %total LH2 volume for ascent/descent + transfer (m^3)
        [A_tank_LLV_LH2,l_tank_LLV_LH2]=tank_area(V_LLV_LH2,d_LLV); %LH2 tank surface area (m^2)
        m_insulation_LLV_LH2=2.88*A_tank_LLV_LH2; %LH2 tank insulation mass (kg)   
        %LLV total wet mass
        %m_LLV_prop=m_LLV_prop_desc+m_LLV_prop_asc+m_bus_prop; %LLV total propellant mass (kg)
        m_LLV_prop=m_LLV_LO2_desc+m_LLV_LO2_asc+m_LLV_LH2_desc+m_LLV_LH2_asc+m_bus_prop; %LLV total propellant mass (kg)
        m_LLV_wet=m_LLV_dry+m_LLV_prop; %LLV wet mass (kg) 
    %fairings
        A_fore=pi*d_LLV*(d_LLV/2); %foreward fairing area (m^2)
        A_intertank=pi*d_LLV*(d_LLV); %intertank fairing area (m^2)
        A_aft=pi*d_LLV*(d_LLV/2); %aft fairing area (m^2)
        A_LLV_fairing=A_fore+A_intertank+A_aft; %total LLV fairing surface area (m^2)
        m_LLV_fairing=4.95*(A_LLV_fairing)^1.15; %mass of LLV fairing (kg)
    %avionics
        m_LLV_avionics=10*m_LLV_wet^0.361; %mass of avionics (kg)
    %wiring
        l_LLV=l_tank_LLV_LO2+l_tank_LLV_LH2; %LLV length (m)
        m_LLV_wiring=1.058*sqrt(m_LLV_wet)*l_LLV^0.25; %mass of wiring (kg)
    %thrust structure
        m_LLV_thrust=n_engine*(2.55e-04)*T; %thrust structure mass (kg)
    
   
    m_LLV_tank_LO2=m_LLV_tank_LO2_desc+m_LLV_tank_LO2_asc+m_bus_tank_LO2+m_insulation_LLV_LO2; %total LLV LO2 tank mass (kg)
    m_LLV_tank_LH2=m_LLV_tank_LH2_desc+m_LLV_tank_LH2_asc+m_bus_tank_LH2+m_insulation_LLV_LH2; %total LLV LH2 tank mass (kg)
    m_LLV_tank=m_LLV_tank_LO2+m_LLV_tank_LH2; %total LLV propellant tank mass (kg)
    dm_tank_LO2=m_LLV_tank_LO2-m_LLV_tank_LO2_estimate; %LO2 tank mass change (kg)
    dm_tank_LH2=m_LLV_tank_LH2-m_LLV_tank_LH2_estimate; %LH2 tank mass change (kg)
    m_LLV_tank_LO2_estimate=m_LLV_tank_LO2; %new calculated mass becomes new mass estimate (kg)
    m_LLV_tank_LH2_estimate=m_LLV_tank_LH2; %new calculated mass becomes new mass estimate (kg)
    m_LLV_dry=m_LLV_crew+m_LLV_engine+m_LLV_COPV+m_legs+m_LLV_tank+m_LLV_fairing+m_LLV_avionics+m_LLV_wiring+m_LLV_thrust; %LLV dry mass (kg)
    count=count+1;
end
disp('Lunar Landing Vehicle')
m_LLV_dry %dry mass of LLV (kg)
m_LLV_prop %LLV propellant mass (kg)
m_LLV_wet %LLV wet mass (kg)
% m_LLV_prop_LO2 %LO2 propellant mass (kg)
% m_LLV_prop_LH2 %LH2 propellant mass (kg)

% LLV initial launch
m_LLV_landing=(m_LLV_dry)*exp((dv_LLV_desc)/(Isp*g0)); %LLV total mass before landing (kg)
m_LLV_prop_landing=m_LLV_landing-m_LLV_dry; %LLV initial descent propellant mass (kg)
m_LLV_launch=(m_LLV_landing)*exp((dv_TLI)/(Isp*g0)); %LLV launch mass (kg)
m_LLV_launch_prop=m_LLV_launch-m_LLV_landing; %LLV ascent propellant mass (kg)
n_tanker_launches=ceil((m_LLV_prop_landing+m_LLV_launch_prop)/m_tanker_prop) %number of tanker flights needed
n_tanker_launches_roundtrip=ceil(m_LLV_prop/m_tanker_prop) %number of launches, full load (kg)
%LLV launch thrust check
F_liftoff=1000*m_LLV_wet*g0_moon %force needed to lift LLV off lunar surface (N)
T_LLV=n_engine*T %LLV liftoff thrust (N)


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