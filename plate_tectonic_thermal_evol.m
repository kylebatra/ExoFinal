function dTdt = plate_tectonic_thermal_evol(t,Y,Q0_forward,mu_ref,Ts,Rp,Rc,Ev)

dTdt=zeros(1,1);
%kyle batra
% Same parameters as listed in "single_model_classic_thermal_evol.m"
% Numbers should be the same between the two files
Rg=8.314;
rho=4000.0;
g=9.8;
d=2890000.0;
alpha=3e-5;
kappa=1e-6;
Cp=1250; % J kg^-1 K^-1
k=5; % W m^-1 K^-1
C1=0.1; % heat flow scaling law constant
beta=1/3;
mun=mu_ref/exp(Ev/(Rg*1650)); 

V_man=(4/3)*pi*(Rp^3-Rc^3);
A_surf=4*pi*Rp^2;

tau_rad2=2.94e9*3600*24*365; % Radioactive decay constant, in seconds 

% Calculate thickness of the lithosphere for plate tectonics
% Define rayleigh number & theta
Rai=(rho*g*alpha*d^3*(Y(1)-Ts))/(kappa*mun*exp(Ev/(Rg*Y(1))));
%theta=(Ev*(Y(1)-Ts))/(Rg*Y(1)^2);
delta=(d/C1)*Rai^(-beta);

% Y(1) is mantle potential temperature; this equation is the right hand
% side of the planetary energy balance equation, and is equal to dT/dt 
dTdt(1)=Q0_forward*exp(-t/tau_rad2)*(1/(V_man*rho*Cp))-(A_surf/(rho*Cp*V_man))*(k*(Y(1)-Ts)/delta);


