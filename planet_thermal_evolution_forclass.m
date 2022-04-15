close all
clear all
%kyle batra, last modifiged 4/15/22
% Units are all m, kg, s 
Ev=300000.0; % Activation energy for viscosity (controls how strongly viscosity depends on temperature)
Rp=1*6378100.0; % Radius of whole planet 
Rc=1*3488.1*1000; % Core radius
Ts=273.0; % Surface temperature (fixed in time)

% Reference viscosity in pascal*seconds; this is the average 
%mantle viscosity at today's temperature; one of the more important variables to change
mu_ref=1e21; % Reference viscosity is read in to solve the thermal evolution equation

V_man=(4/3)*pi*(Rp^3-Rc^3);
Q0_forward_Earth=100*1e12; % Initial radioactive heating power in TW for Earth
V_man_Earth=9.06e20; % Volume of Earth's mantle in m^3
Q_rad=Q0_forward_Earth/V_man_Earth;
Q0_forward=Q_rad*V_man;
tau_rad2=2.94e9*3600*24*365; % In seconds

t_end=1e10*3600*24*365; % time when we stop the model; 4.5e9 is in years, then converted into seconds
Ti_init=2000; % Initial mantle temperature

% These you can leave as is; they are the tolerances for the ode solver
options = odeset('RelTol',1e-13,'AbsTol',1e-10);
% This calls matlab's ode solver to solve for the time evolution of
% temperature
[T,Y]=ode15s(@(t,y) stagnant_thermal_evol(t,y,Q0_forward,mu_ref,Ts,Rp,Rc,Ev),...
    [0 t_end],Ti_init,options);

time(1:length(T))=T(1:length(T))/(1e9*3600*24*365); % convert time from seconds to billions of years
Ti(1:length(T))=Y(1:length(T),1); % Take the solution for mantle temperature as a function of time, put in vector Ti

[N,M]=ode15s(@(t,y) plate_tectonic_thermal_evol(t,y,Q0_forward,mu_ref,Ts,Rp,Rc,Ev),...
    [0 t_end],Ti_init,options);

time2(1:length(N))=N(1:length(N))/(1e9*3600*24*365); % convert time from seconds to billions of years
Ti2(1:length(N))=M(1:length(N),1); % Take the solution for mantle temperature as a function of time, put in vector Ti
tectlim=1120.661;
% Plots mantle potential temperature as a function of time
figure(1)
plot(time,Ti,time2,Ti2)
xlabel('Time [Gyrs]')
ylabel('Mantle Temperature [K]');
legend('Stagnant Lid','Tectonic Plate','location','best')
%% calculate heat loss and heat production over time

times(1:length(T))=T(1:length(T));
Q_rad=Q0_forward*exp(-times/tau_rad2); 

times2(1:length(N))=N(1:length(N));
Q_rad2=Q0_forward*exp(-times/tau_rad2);

Ev=300000.0;
Rg=8.314;
rho=4000.0;
g=9.8;
d=2890000.0;
alpha=3e-5;
kappa=1e-6;
Rp=6378100.0;
Cp=1250; % J kg^-1 K^-1
k=5; % W m^-1 K^-1
C1s=0.5; % heat flow scaling law constant staglid
C1=0.1; % heat flow scaling law constant
beta=1/3;
mun=mu_ref/exp(Ev/(Rg*1650)); 

V_man=(4/3)*pi*(Rp^3-Rc^3);
A_surf=4*pi*Rp^2;

% Calculate thickness of the lithosphere for plate tectonics
% Define rayleigh number & theta
Rai=(rho*g*alpha*d^3*(Ti-Ts))./(kappa*mun*exp(Ev./(Rg*Ti)));
Rai2=(rho*g*alpha*d^3*(Ti2-Ts))./(kappa*mun*exp(Ev./(Rg*Ti2)));
theta=(Ev.*(Ti-Ts))./(Rg*Ti.^2); %y(1)has to be ti
delta=(d/C1s)*theta.^(4/3).*Rai.^(-beta);
delta2=(d/C1)*Rai2.^(-beta);
heat_loss=A_surf*(k*(Ti-Ts)./delta);
heat_loss2=A_surf*(k*(Ti2-Ts)./delta2);

pressure=(g.*delta*rho)./(10^(9));
MeltTemp=-5.104.*pressure.^(2)+132.899.*pressure+1120.661+273.15;
MeltTempTP=((1120.661+273.15).*delta2)./delta2;

figure(2)
plot(time,delta,time2,delta2)
xlabel('Time [Gyrs]') 
ylabel('Delta Thickness [M]')
legend('Stagnant Lid Thickness','Plate Tectonic Crust Thickness','location','best')


figure(3)
plot(time,heat_loss/1e12,time2,heat_loss2/1e12,time,Q_rad/1e12)
xlabel('Time [Gyrs]') 
ylabel('Heat production/Heat loss [TW]')
legend('Mantle heat loss (Stagnant Lid)','Mantle heat loss (Plate Tectonics)','Radioactive heat production','location','best')

figure(4)
plot(time,Ti,time2,Ti2,time,MeltTemp,time2,MeltTempTP)
xlabel('Time [Gyrs]')
ylabel('Mantle Temperature [K]');
legend('Stagnant Lid Interior Temperature','Tectonic Plate Interior Temperature','Stagnant Lid Melting Temp','Tectonic Plate Melting Temp','location','best')