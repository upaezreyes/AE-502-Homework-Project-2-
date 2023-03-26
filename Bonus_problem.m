%AE 502 Bonus problem
% Gauss Planetary Equatons of equinoctial elements 
% using the method of Lagrange/Poisson brackets 

function Bonus_problem

clc
clear variables 

pi = 3.1416; % pi
mu = 3.986e5; % km^3/s^2, gravitational constant
Re = 6378; % km, Earth's radius
J2 = 0.00108; % J2 perturbation value

% Given Equinoctial Elements: 
a0 = 7000; % km, initial, semi-major axis
e0 = 0.05; % initial, eccentricity 
%i0 = 0 ; % radians, initial, inclination angle
h0 = 0;  % h , initial value
p0 = 0;  % p , iniitial value
q0 = 0;  % q , iniitial value
k0 = e0; % k , iniitial value
%p0 = a0*(1-e0.^2); % km, initial semi-latus rectum 
omega0 = atan(p0./q0); % radians, right ascension of the node
%i0 = 2.*atan(sqrt((p0.^2) + (q0.^2))); % radians, inclination angle
w0 = atan(h0./k0) - atan(p0./q0);
f0 = 0; 
L0 = w0 + omega0 + f0; 
l0 = 0; 
 

T0 = 2.*pi.*sqrt((a0.^3)./mu); % sec, initial period

% store initial equinoctial elements: 
coe0 = [a0 h0 p0 q0 k0 L0 l0]; 

% Use ODE45 to integrate the Gauss variational equations (Equations
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

% Assign the time histories mnemonic variable names:
a = y(:,1); 
h = y(:,2); 
p = y(:,3); 
q = y(:,4); 
k = y(:,5); 
L = y(:,6); 
l = y(:,7); 

% for gh=1:length(y)
%     fprintf('p = %f \n',p)
% end
e = sqrt((h.^2) + (k.^2)); % eccentricity
%a = p./(1 - e.^2); % km, semi-major axis
omega = atan(p./q); % radians, right ascension of the node
i = 2.*atan(sqrt((p.^2) + (q.^2))); % radians, inclination angle
w = atan(h./k) - atan(p./q); % radians, argument of perigee
%w = asin(h./e) - omega;
f = L - atan(h./k); % radians, anomaly angle

%disp(a)
disp(e)

%Plot the time histories of the osculating elements:
figure
subplot(2,1,1)
plot(t/3600,omega*(180./pi))
title('Right Ascension (degrees)')
xlabel('hours')
ylabel('degrees')
xlim([0 24])
grid on
subplot(2,1,2)
plot(t/3600,w*(180./pi))
title('Argument of Perigee (degrees)')
xlabel('hours')
ylabel('degrees')
xlim([0 24])
grid on


figure
subplot(2,1,1)
plot(t/3600,e)
title('Eccentricity')
xlabel('hours')
ylabel('e')
grid on
subplot(2,1,2)
plot(t/3600,a)
title('Semi-major Axis (km)')
xlabel('hours')
ylabel('km')
grid on

figure
subplot(2,1,1)
plot(t/3600,i*(180./pi))
title('Inclination')
xlabel('hours')
ylabel('degrees')
grid on
subplot(2,1,2)
plot(t/3600,L*(180./pi))
title('True Longitude (degrees)')
xlabel('hours')
ylabel('degrees')
grid on

%=========================================================================
% Subfunction: 
function dfdt = rates(~,f)

% equinoctial elements at time t: 
a = f(1); 
h = f(2); 
p = f(3); 
q = f(4); 
k = f(5); 
L = f(6); 
l = f(7); 

e = sqrt((h.^2) + (k^2)); % eccentricity
%a = p./(1 - e.^2); % km, semi-major axis
omega = atan(p./q); % radians, right ascension of the node
i = 2.*atan(sqrt((p.^2) + (q.^2))); % radians, inclination angle
w = atan(h./k) - atan(p./q); % radians, argument of perigee
%w = asin(h./e) - omega;
f = L - atan(h./k); % radians, true anomaly 
u = w + f; % radians, argument of latitude 

% constants: 
n = sqrt(mu./(a.^3)); 
b = a.*sqrt(1 - (h.^2) - (k.^2)); 
hbar = n.*a.*b; 
pbar_r = 1 + h.*sin(L) + k.*cos(L); 
r_hbar = hbar./(mu.*(1 + h.*sin(L) + k.*cos(L)));  
J = (hbar.^2)./mu; 
r = J./(1 + h.*sin(L) + k.*cos(L)); 

%perturbation accelerations: 
% radial direction: 
adr = -(3./2).*J2.*mu.*((Re.^2)./r.^4).*(1-3.*(sin(i).^2).*(sin(u).^2)); 

% ortogonal/transverse direction: 
adt = -(3./2).*J2.*mu.*((Re.^2)./r.^4).*(sin(i).^2).*sin(2.*u); 

% angular momentum direction: 
adh = -(3./2).*J2.*mu.*((Re.^2)./r.^4).*sin(2.*i).*sin(u); 

%===
% BB = sqrt(1 - h.^2 - k.^2); 
% phi = 1 + h.*sin(L) + k.*cos(L); 

% equinoctial element rates: 
% da/dt: 
adot = ((2.*a.^2)./hbar).*((k.*sin(L) - h.*cos(L)).*adr + pbar_r.*adt); 
%adot = (2./BB).*sqrt((a.^3)./mu).*((k.*sin(L) - h.*cos(L)).*adr + phi.*adt); 

% dh/dt: 
hdot = r_hbar.*(-(pbar_r).*cos(L).*adr + (h + (1 + pbar_r).*sin(L)).*adt ...
        - k.*(p.*cos(L) - q.*sin(L)).*adh); 
%hdot = BB.*sqrt(a./mu).*(-adr.*cos(L) + (((h+sin(L))./phi)+sin(L)).*adt - k.*((p.*cos(L) - q.*sin(L))./phi).*adh);

% dk/dt: 
kdot = r_hbar.*(pbar_r.*sin(L).*adr + (k + (1 + pbar_r).*cos(L)).*adt ...
        + h.*(p.*cos(L) - q.*sin(L)).*adh);
%kdot = BB.*sqrt(a./mu).*(adr.*cos(L) + (((k+cos(L))./phi)+cos(L)).*adt + h.*((p.*cos(L) - q.*sin(L))./phi).*adh);

% dp/dt: 
pdot = (1./2).*r_hbar.*(1 + p.^2 + q.^2).*sin(L).*adh; 
%pdot = (BB./2).*sqrt(a./mu).*(1 + p.^2 + q.^2).*(sin(L)./phi).*adh; 

% dq/dt: 
qdot = (1./2).*r_hbar.*(1 + p.^2 + q.^2).*cos(L).*adh; 
%qdot = (BB./2).*sqrt(a./mu).*(1 + p.^2 + q.^2).*(cos(L)./phi).*adh; 

% dL/dt: 
% B = 1 + h.*sin(L) + k.*cos(L); 
% C = sqrt(1 - h.^2 - k.^2); 

Ldot = (((mu.^2).*(1 + h.*sin(L) + k.*cos(L)).^2)./hbar) ...
        + (hbar./mu).*((q.*sin(L) - p.*cos(L))./(1 + h.*sin(L) + k.*cos(L))).*adh; 
% Ldot = sqrt(mu./a).*((phi.^2)./BB.^3) ... 
%         - sqrt((a.^3)./mu).*(BB./phi).*adh.*(p.*cos(L) - q.*sin(L)); 
% 

%dl/dt: 
l1 = (a./(a+b)).*(pbar_r).*(h.*sin(L)+k.*cos(L)) + ((2.*b)./a); 
l2 = (a./(a+b)).*(1+pbar_r).*(h.*cos(L)-k.*sin(L)); 
l3 = p.*cos(L) - q.*sin(L); 
ldot = n - r_hbar.*(l1.*adr + l2.*adt + l3.*adh); 

% % constants: 
% W = 1 + fa.*cos(L) + g.*sin(L); 
% s2 = 1 + (h.^2) + (k.^2); 
% r = p./W; % km, distance 
% 
% % components of the perturbation acceleration: 
% % radial component ft: 
% fr = (3./2).*mu.*J2.*((Re.^2)./r.^4).*(3.*(sin(i).^2).*(sin(u).^2) - 1); 
% % transverse component ft: 
% ft = -3.*mu.*J2.*((Re.^2)./r.^4).*(sin(i).^2).*sin(u).*cos(u); 
% % normal component: 
% fn = -3.*mu.*J2.*((Re.^2)./r.^4).*sin(i).*cos(i).*sin(u); 
% 
% % equinoctial elements: 
% % dp/dt: 
% pdot = ((2.*p)./W).*sqrt(p./mu).*ft; 
% 
% % df/dt: 
% fdot = sqrt(p./mu).*(fr.*sin(L)+((W+1).*cos(L)+fa).*(ft./W)  ...
%        - (h.*sin(L)-k.*cos(L)).*(g./W).*fn); 
% 
% % dg/dt: 
% gdot = sqrt(p./mu).*(-fr.*cos(L)+((W+1).*sin(L)+g).*(ft./W) ...
%        + (h.*sin(L)-k.*cos(L)).*(fa./W).*fn); 
% 
% % dh/dt: 
% hdot = sqrt(p./mu).*((s2.*fn)./(2.*W)).*cos(L); 
% 
% % dk/dt: 
% kdot = sqrt(p./mu).*((s2.*fn)./(2.*W)).*sin(L); 
% 
% % dL/dt: 
% Ldot = sqrt(mu.*p).*((W./p).^2) ...
%        + (1./W).*sqrt(p./mu).*(h.*sin(L)-k.*cos(L)).*fn; 

% pass rates back to ODE45 in the array dfdt: 
dfdt = [adot hdot pdot qdot kdot Ldot ldot]'; 
end 

end 




