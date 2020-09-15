%John Furumo
%ENAE 791 section 0101
%Lunar Mission Project
%date 04/02/2020

clear;
clf;
clc;
close all;
format shortg;
set(0,'DefaultFigureWindowStyle','docked');

%% CONSTANTS
G=6.67430e-20; %km^3/kg*s^2

%Earth
m_earth=5.97219e+24; %kg
mu_earth=398600.4415; %km^3/s^2
R_earth=6378; %km

%Moon
m_moon=7.349e+22; %kg
mu_moon=4902.800066; %km^3/s^2 
R_moon=1737.4; %km


%% Orbits
%Moon orbit
a_moon=384400; %km
e_moon=0.05490;
r_p_moon=a_moon*(1-e_moon); %km
r_a_moon=a_moon*(1+e_moon); %km
i_moon=5.145; %deg
P_moon=2*pi*sqrt(a_moon^3/mu_earth); %s
%Moon sphere of influence
r_SOI_moon=a_moon*(mu_moon/mu_earth)^(2/5); %km

%LEO parking orbit
a_LEO=10000; %km
e_LEO=0;
i_LEO=28.5; %deg
P_LEO=2*pi*sqrt(a_LEO^3/mu_earth); %s
v_LEO=sqrt(mu_earth/a_LEO); %km/s

%TLI
r_p_TLI=a_LEO; %km
r_a_TLI=r_p_moon-r_SOI_moon; %km
a_TLI=(r_p_TLI+r_a_TLI)/2; %km
e_TLI=r_a_TLI/a_TLI-1;
i_TLI=i_moon; %deg
P_TLI=2*pi*sqrt(a_TLI^3/mu_earth); %s
v_TLI_p=sqrt(2*mu_earth/r_p_TLI-mu_earth/a_TLI); %km/s
v_TLI_a=sqrt(2*mu_earth/r_a_TLI-mu_earth/a_TLI); %km/s

%LLO parking orbit
h_LLO=100; %km
a_LLO=h_LLO+R_moon; %km
e_LLO=0;
i_LLO=90; %deg
P_LLO=2*pi*sqrt(a_LLO^3/mu_moon); %s
v_LLO=sqrt(mu_moon/a_LLO); %km/s



%% Time domains

%Launch
t_init_launch=0; %s
t_res_launch=60; %s
t_term_launch=600; %s
t_span_launch=[t_init_launch:t_res_launch:t_term_launch]; %s

%Earth parking orbit
t_init_park=t_term_launch; %s
t_res_park=60; %s
t_term_park=3600; %s
t_span_park=[t_init_park:t_res_park:t_term_park]; %s

%TLI - transfer trajectory
t_init_TLI=t_term_park; %s
t_res_TLI=60; %s
t_term_TLI=3600; %s
t_span_TLI=[t_init_TLI:t_res_TLI:t_term_TLI]; %s

%LOI - lunar parking orbit
t_init_LOI=t_term_TLI; %s
t_res_LOI=60; %s
t_term_LOI=3600; %s
t_span_LOI=[t_init_LOI:t_res_LOI:t_term_LOI]; %s


%Lunar orbit plane change

%Descent

%Ascent

%Lunar orbit

%TEI

%Entry





delta_v_1=sqrt(v_LEO^2+v_TLI_p^2-2*v_LEO*v_TLI_p*cosd(i_LEO-i_moon)) %km/s

delta_v_inc=2*v_LEO*sind((i_LEO-i_moon)/2);

delta_v_HT=v_TLI_p-v_LEO;


delta_v_TEI=sqrt(v_LLO^2+v_TLI_a^2-2*v_LLO*v_TLI_a*cosd(i_LLO-i_moon)) %km/s


%% LOI

%Maneuver 2: Transfer HEO - LLO
ToF=P_TLI/2; %s

i_LLO=i_moon;


v_TLI_a=sqrt(2*mu_earth/r_a_TLI-mu_earth/a_TLI) %km/s

epsilon_LOI=v_TLI_a^2/2-mu_moon/r_SOI_moon %km^2/s^2

% oe_LOI=[r_SOI_moon;0;i_TLI;0;0;180] %OE state vector
% cart_LOI=oe2cart(oe_LOI)

% cart_LOI=rv_TLI(t_TLI(end)/2,:)'

%% Numerical integration orbits (2BP)
%P=27.321582; %days

oe_moon=[a_moon;e_moon;i_moon;0;0;0]; %OE state vector
cart_moon = oe2cart(oe_moon,mu_earth); %cartesian state vector

oe_park=[a_park;e_park;i_park;0;0;0]; %OE state vector
cart_park=oe2cart(oe_park,mu_earth);

oe_TLI_0=[a_TLI;e_TLI;i_moon;0;180;0] %OE state vector
cart_TLI_0=oe2cart(oe_TLI_0,mu_earth);


%parking orbit
P=28; %days
P_park=P*86400; %s
t_span_park=[0:60:P_park]; %s
tolerance=1e-009;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
%moon orbit
[t_park,rv_moon1]=ode45(@integrate_2BP,t_span_park,cart_moon,options,mu_earth);
%mission parking orbit
[t_park,rv_park]=ode45(@integrate_2BP,t_span_park,cart_park,options,mu_earth);

%TLI
TLI_end=P_park+P_TLI/2;
t_span_TLI=[P_park:60:TLI_end];
v_final = rotate_vector(rv_park(end,4:6)',deg2rad(-(i_park-i_moon)),rv_park(end,1:3)');

cart_TLI_0_B=[rv_park(end,1:3)';v_TLI_p*v_final/norm(v_final)]
oe_TLI_0=cart2oe(cart_TLI_0_B,mu_earth)
%moon orbit
[t_TLI,rv_moon2]=ode45(@integrate_2BP,t_span_park,rv_moon1(end,:),options,mu_earth);
%trans-lunar injection
[t_TLI,rv_TLI]=ode45(@integrate_2BP,t_span_park,cart_TLI_0,options,mu_earth);

t_total=t_park+t_TLI; %s

%LOI
cart_LOI=[[r_p_moon;0;0]-rv_TLI(20161,1:3)';rv_TLI(20161,4:6)']
t_span_LOI=[TLI_end:60:TLI_end+86400];
[t_LOI,rv_LOI]=ode45(@integrate_2BP,t_span_LOI,cart_LOI,options,mu_moon);


% %% Animated Plot, Earth - Moon system
% figure('name','Earth - Moon System');
% axis equal;
% rotate3d on
% res=8;
% %view([500000,500000,0000]);
% for i=1:12:length(t_park)
%     %plot Earth
% %     loc_earth=[0,0,0];
% %     p_earth=plot_sphere(loc_earth,radius_earth,res,'Earth','b');
%     
%     plot3(0,0,0,'bo')
%     hold on;
%     %plot Moon orbit
%     p_moon_orbit = plot3(rv_moon1(:,1),rv_moon1(:,2),rv_moon1(:,3),'g-','linewidth',1,'DisplayName','Lunar Orbit');
%     %plot Moon
%     %loc_moon=[rv_moon(i,1),rv_moon(i,2),rv_moon(i,3)];
%     %p_moon=plot_sphere(loc_moon,r_moon,res,'Moon','k');
%     %plot Moon SOI
%     %p_moon_SOI=plot_sphere(loc_moon,r_SOI_moon,16,'Moon Sphere of Influence','r');
%     p_moon=plot3(rv_moon1(i,1),rv_moon1(i,2),rv_moon1(i,3),'ko');
%     %plot mission trajectory
%     p_park=plot3(rv_park(i,1),rv_park(i,2),rv_park(i,3),'rs');
%     plot3(rv_park(:,1),rv_park(:,2),rv_park(:,3),'r-');
%     pause(0.1)
%     hold off
%     axis equal;
% %     title('Earth - Moon System');
%     xlabel('X Position (km)');
%     ylabel('Y Position (km)');
%     zlabel('Z Position (km)');
%     %legend([p_earth p_moon_orbit p_moon p_park])
%     %drawnow
% end
% 
% for i=1:12:length(t_TLI)
%    
%     plot3(0,0,0,'bo')
%     hold on;
%     %plot Moon orbit
%     p_moon_orbit = plot3(rv_moon1(:,1),rv_moon1(:,2),rv_moon1(:,3),'g-','linewidth',1,'DisplayName','Lunar Orbit');
%     %plot Moon
%     p_moon=plot3(rv_moon2(i,1),rv_moon2(i,2),rv_moon2(i,3),'ko');
%     %plot mission trajectory
%     p_TLI=plot3(rv_TLI(i,1),rv_TLI(i,2),rv_TLI(i,3),'rs');
%     plot3(rv_TLI(:,1),rv_TLI(:,2),rv_TLI(:,3),'r-');
%     pause(0.1)
%     hold off
%     axis equal;
% %     title('Earth - Moon System');
%     xlabel('X Position (km)');
%     ylabel('Y Position (km)');
%     zlabel('Z Position (km)');
%     %legend([p_earth p_moon_orbit p_moon p_park])
%     %drawnow
% end

%% Static Plot, Earth-Moon System
% figure('name','Parking Orbit')
% axis equal;
% rotate3d on
% plot3(0,0,0,'bo')
% hold on;
% %plot Moon orbit
% p_moon_orbit = plot3(rv_moon1(:,1),rv_moon1(:,2),rv_moon1(:,3),'g-','linewidth',1,'DisplayName','Lunar Orbit');
% %plot Moon
% res=8;
% loc_moon=[rv_moon1(end,1),rv_moon1(end,2),rv_moon1(end,3)];
% p_moon=plot_sphere(loc_moon,r_moon,res,'Moon','k');
% %plot Moon SOI
% p_moon_SOI=plot_sphere(loc_moon,r_SOI_moon,16,'Moon Sphere of Influence','r');
% p_moon=plot3(rv_moon1(end,1),rv_moon1(end,2),rv_moon1(end,3),'ko');
% %plot mission trajectory
% p_park=plot3(rv_park(end,1),rv_park(end,2),rv_park(end,3),'rs');
% plot3(rv_park(:,1),rv_park(:,2),rv_park(:,3),'r-');
% hold off
% axis equal;
% title('Earth - Moon System');
% xlabel('X Position (km)');
% ylabel('Y Position (km)');
% zlabel('Z Position (km)');


figure('name','TLI')
res=8;
axis equal;
rotate3d on
plot3(0,0,0,'bo')
hold on;
%plot Moon orbit
p_moon_orbit = plot3(rv_moon1(:,1),rv_moon1(:,2),rv_moon1(:,3),'g-','linewidth',1,'DisplayName','Lunar Orbit');
[~,~] = plot_vector_3D(rv_moon1(1,1:3),[0,0,0],'b','-','x')

%plot Moon
loc_moon=[rv_moon1(1,1),rv_moon1(1,2),rv_moon1(1,3)];
p_moon=plot_sphere(loc_moon,r_moon,res,'Moon','k');
%plot Moon SOI
p_moon_SOI=plot_sphere(loc_moon,r_SOI_moon,16,'Moon Sphere of Influence','r');
p_moon=plot3(rv_moon1(1,1),rv_moon1(1,2),rv_moon1(1,3),'ko');
%plot mission trajectory
plot3(rv_park(:,1),rv_park(:,2),rv_park(:,3),'r-');
p_TLI=plot3(rv_TLI(1,1),rv_TLI(1,2),rv_TLI(1,3),'rs');
plot3(rv_TLI(:,1),rv_TLI(:,2),rv_TLI(:,3),'r-');
hold off
axis equal;
title('Earth - Moon System');
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');

figure('name','LOI')
res=8;
rotate3d on
%plot Moon
plot3(0,0,0,'ko')
hold on;
loc_moon=[0,0,0];
p_moon=plot_sphere(loc_moon,r_moon,res,'Moon','k');
%plot Moon SOI
p_moon_SOI=plot_sphere(loc_moon,r_SOI_moon,16,'Moon Sphere of Influence','r');
%plot mission trajectory
p_LOI=plot3(rv_LOI(1,1),rv_LOI(1,2),rv_LOI(1,3),'rs');
plot3(rv_LOI(:,1),rv_LOI(:,2),rv_LOI(:,3),'r-');
hold off
axis equal;
title('Moon System');
xlabel('X Position (km)');
ylabel('Y Position (km)');
zlabel('Z Position (km)');




% figure('name','velocity vector check')
% 
% rotate3d on
% plot3(0,0,0,'bo')
% hold on;
% plot3(rv_park(:,1),rv_park(:,2),rv_park(:,3),'r-');
% plot3(rv_TLI(:,1),rv_TLI(:,2),rv_TLI(:,3),'c-');
% %p_moon_orbit = plot3(rv_moon2(:,1),rv_moon2(:,2),rv_moon2(:,3),'c-','linewidth',1,'DisplayName','Lunar Orbit');
% [~,~] = plot_vector_3D(rv_park(end,1:3),[0,0,0],'b','-','x')
% [~,~] = plot_vector_3D(1000*rv_park(end,4:6),rv_park(end,1:3),'g','-','x')
% [~,~] = plot_vector_3D(1000*cart_TLI_0(4:6),rv_park(end,1:3),'g','-','x')
% hold off;
% axis equal;
% view([rv_park(end,1),rv_park(end,2),rv_park(end,3)])
% xlabel('X Position (km)');
% ylabel('Y Position (km)');
% zlabel('Z Position (km)');

[dv1,dv2,ToF] = Hohmann_Transfer_ind(r_SOI_moon,100+r_moon,mu_moon)
dv_LOI_lower=dv1+dv2
v_LLO=sqrt(mu_moon/(r_moon+100))
dv_LOI_plane_change=2*v_LLO*sind((90)/2)
dv_LOI_comb=sqrt(v_LLO^2+v_TLI_a^2-2*v_LLO*v_TLI_a*cosd(90))

%"Hohmann Transfer" launch trajectory
r1=r_earth;
r2=r_park;
mu=mu_earth;
v_equator=2*pi*r_earth/86400 %%km/s
v_c1=v_equator; %km/s
%trangential velocity of transfer ellipse at initial apse
a=(r1+r2)/2; %km
v_t1=sqrt(2*mu/r1-mu/a); %km/s
%first burn delta-v
dv1=abs(v_t1-v_c1); %km/s
%tangential velocity of final orbit
v_c2=sqrt(mu/r2); %km/s
%tangential velocity of transfer ellipse at final apse
v_t2=sqrt(2*mu/r2-mu/a); %km/s
%final burn delta-v
dv2=abs(v_t2-v_c2); %km/s
%total delta-v
dv_Hohmann=dv1+dv2 %km/s

