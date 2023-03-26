function Example_12_6

clc
clear variables 

pi = 3.1416; % pi
mu = 3.986e5; % km^3/s^2, gravitational constant
Re = 6378; % km, Earth's radius 
J2 = 0.00108; % J2 pertubation 

% Given Orbital Elemens: 
% a0 = 26600; % km, initial semi-major axis 
% i0 = .4887; % degrees, initial inclination angle 
% e0 = 0.17136; % initial eccentricity 
% a0 = 6678/(1-e0);
% w0 = 30.*(pi./180); % radians, initial argument of perigee angle
% omega0 = 45.*(pi./180); % radians, initial right ascension of the node
% M0 = 40.*(pi./180); % radians, initial man anomaly 

%a0 = 26600; % km, initial semi-major axis 
i0 = 28.*(pi./180); % degrees, initial inclination angle 
rp0 = 6678; 
ra0 = 9940; 
e0 = (ra0-rp0)./(ra0+rp0); % initial eccentricity 
a0 = rp0/(1-e0);
w0 = 30.*(pi./180); % radians, initial argument of perigee angle
omega0 = 45.*(pi./180); % radians, initial right ascension of the node
M0 = 40.*(pi./180); % radians, initial man anomaly 

%rp0 = a0.*(1-e0); % km, initial radius of perigee
%ra0 = a0.*(1+e0); % km, initial radius of apogee 
h0 = sqrt(mu.*a0.*(1-e0.^2)); % km^2/s, initial angular momentum 
T0 = 2.*pi.*sqrt((a0.^3)./mu); % sec, initial period 

% store initial orbital elements: 
coe0 = [h0 e0 omega0 i0 w0 M0]; 


%...Use ODE45 to integrate the Gauss variational equations (Equations
% 12.89) from t0 to tf:
t0 = 0; % sec, initial time
days = 2; % days of time evolution
tf = days*24*3600; % sec, final time
nout = 5000; %number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout); % time steps
options = odeset(...
 'reltol', 1.e-8, ... 
 'abstol', 1.e-8, ...
'initialstep', T0/1000);

y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);

%...Assign the time histories mnemonic variable names:
h = y(:,1);
e = y(:,2);
omega = y(:,3);
i = y(:,4);
w = y(:,5);
M = y(:,6);

a = (h.^2)./(mu.*(1-e.^2));

%...Plot the time histories of the osculating elements:
figure
subplot(6,1,1)
plot(t/3600,omega*(180/pi))
title('Right Ascension (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(6,1,2)
plot(t/3600,w*(180/pi))
title('Argument of Perigee (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(6,1,3)
plot(t/3600,h - h0)
title('Angular Momentum (km^2/s)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(6,1,4)
plot(t/3600,e)
title('Eccentricity')
xlabel('hours')
grid on
grid minor
axis tight

subplot(6,1,5)
plot(t/3600,i*(180/pi))
title('Inclination (degrees)')
xlabel('hours')
grid on
grid minor
axis tight

subplot(6,1,6)
plot(t/3600, a)
title('Semi-major Axis (km)')
xlabel('hours')
ylabel('km')
grid on

% subfunction:
%====================================================================
function dfdt = rates(t,f)
%w===================================================================
%
% This function calculates the time rates of the orbital elements
% from Gauss’s variational equations (Equations 12.89).
%-------------------------------------------------------------------------
%...The orbital elements at time t:
h = f(1);
e = f(2);
omega = f(3);
i = f(4);
w = f(5);
M = f(6);

r = h^2/mu/(1 + e*cos(M)); % the radius
u = w + M; % argument of latitude

% orbital element rates at time t (Equations 12.89):
hdot = -3/2*J2*mu*Re^2/r^3*sin(i)^2*sin(2*u);

% edot = ...
% 3/2*J2*mu*Re^2/h/r^3*(h^2/mu/r ...
% *(sin(u)*sin(i)^2*(3*sin(M)*sin(u) - 2*cos(M)*cos(u)) - sin(M)) ...
% -sin(i)^2*sin(2*u)*(e + cos(M)));

edot = 3/2*J2*mu*Re^2/h/r^3 ...
*(h^2/mu/r*sin(M)*(3*sin(i)^2*sin(u)^2 - 1) ...
-sin(2*u)*sin(i)^2*((2+e*cos(M))*cos(M)+e));

Mdot = h/r^2 + 3/2*J2*mu*Re^2/e/h/r^3 ...
*(h^2/mu/r*cos(M)*(3*sin(i)^2*sin(u)^2 - 1) ...
+ sin(2*u)*sin(i)^2*sin(M)*(h^2/mu/r + 1));

omegadot = -3*J2*mu*Re^2/h/r^3*sin(u)^2*cos(i);

idot = -3/4*J2*mu*Re^2/h/r^3*sin(2*u)*sin(2*i);

wdot = 3/2*J2*mu*Re^2/e/h/r^3 ...
*(-h^2/mu/r*cos(M)*(3*sin(i)^2*sin(u)^2 - 1) ...
- sin(2*u)*sin(i)^2*sin(M)*(2 + e*cos(M)) ...
+ 2*e*cos(i)^2*sin(u)^2);

% pass these rates back to ODE45 in the array dfdt:
dfdt = [hdot edot omegadot idot wdot Mdot]';
end

end 
