clear;
clc;
format longg;

%{
Main function. Contains code for Q1, Q2, and Q3.
Utilizes functions: oe_to_rv, rv_to_oe, rv_ode
%}


% (Q1)
J2_e = 0.00108;
R_e = 6370; %km
mu_e = 3.986e5; %km^3/s^2
rp_min1 = 600 + R_e; %km

%calculate period (to orbit 3x in one day)
T1 = (1*24*60*60)/3; %s

%calculate semimajor axis
a1 = (mu_e*T1^2/(4*pi^2))^(1/3); %km

%calculate mean motion
n1 = sqrt(mu_e/a1^3); %1/s or rad/s

%choose i s.t. orbit is frozen
i1 = acos(sqrt(1/5)); %rad
i12 = acos(-sqrt(1/5)); %rad

%calculate maximum e
e1 = 1 - (rp_min1)/a1;

%calculate whether orbit is frozen (should = 0)
fr1 = (3/4)*n1*J2_e*(R_e/a1)^2*((5*cos(i1)^2 - 1)/(1-e1^2)^2); %rad/s

%calculate drift rate/nodal precession
dr1 = -(3/2)*n1*J2_e*(R_e/a1)^2*(cos(i1)/(1-e1^2)^2); %rad/s
dr12 = -(3/2)*n1*J2_e*(R_e/a1)^2*(cos(i12)/(1-e1^2)^2); %rad/s

disp('QUESTION 1------------')
disp(['a: ', num2str(a1), ' km'])
disp(['e: ', num2str(e1)])
disp(['i: ', num2str(i1*180/pi), ' deg, ', num2str(i12*180/pi), ' deg'])
disp(['nodal precession: ', num2str(dr1), ' rad/s, ' num2str(dr12), ' rad/s'])



% (Q2)
J2_m = 0.00196;
R_m = 3390; %km
mu_m = 4.282e4; %km^3/s^2
rp_min2 = 400 + R_m; %km

%calculate period (to orbit 1x in one Martian day)
T2 = (24*60*60)+(39*60)+35; %s

%calculate semimajor axis
a2 = (mu_m*T2^2/(4*pi^2))^(1/3); %km

%calculate mean motion
n2 = sqrt(mu_m/a2^3); %1/s or rad/s

%choose i s.t. orbit is frozen
i2 = acos(sqrt(1/5)); %rad
i22 = acos(-sqrt(1/5)); %rad

%calculate maximum e
e2 = 1 - (rp_min2)/a2;

%calculate whether orbit is frozen (should = 0)
fr2 = (3/4)*n2*J2_m*(R_m/a2)^2*((5*cos(i2)^2 - 1)/(1-e2^2)^2); %rad/s

%calculate drift rate/nodal precession
dr2 = -(3/2)*n2*J2_m*(R_m/a2)^2*(cos(i2)/(1-e2^2)^2); %rad/s
dr22 = -(3/2)*n2*J2_m*(R_m/a2)^2*(cos(i22)/(1-e2^2)^2); %rad/s

disp('QUESTION 2------------')
disp(['a: ', num2str(a2), ' km'])
disp(['e: ', num2str(e2)])
disp(['i: ', num2str(i2*180/pi), ' deg, ', num2str(i22*180/pi), ' deg'])
disp(['nodal precession: ', num2str(dr2), ' rad/s, ' num2str(dr22), ' rad/s'])



% (Q3) input data
a = 26600; %km
i = 1.10654; %rad
e = 0.74;
w = 5; %deg
omega = 90; %deg
M0 = 10; %deg

time = 100; %days


% solving using special perturbations equations
%convert to r,v
[r0, v0] = oe_to_rv(a, e, i*180/pi, omega, w, M0, mu_e);
%propagate orbit
duration = time*24*60*60; %s
state0 = [r0 v0]';
tspan = linspace(0, duration, 100000);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,trajectory] = ode45(@rv_ode,tspan,state0,opts,mu_e,J2_e,R_e);
tdays = t/(24*60*60); %days

%plot orbit
xvals = trajectory(:, 1);
yvals = trajectory(:, 2);
zvals = trajectory(:, 3);
figure(1)
hold on
grid on
axis equal
box on
plot3(xvals, yvals, zvals)
plot3(0.,0.,0.,'.','MarkerSize',20) %plot Earth
title('Molniya Orbit about Earth with J2 Perturbation')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(-50,15)
legend('Molniya Orbit','Earth')

%plot orbit another way
figure(2)
hold on
grid on
axis equal
box on
x1 = xvals(1:50000);
y1 = yvals(1:50000);
z1 = zvals(1:50000);
x2 = xvals(50000:100000);
y2 = yvals(50000:100000);
z2 = zvals(50000:100000);
plot3(x1,y1,z1)
plot3(x2,y2,z2)
plot3(0.,0.,0.,'.','MarkerSize',20) %plot Earth
title('Molniya Orbit about Earth with J2 Perturbation')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
view(-50,15)
legend('Orbit (days 1-50)','Orbit (days 50-100)','Earth')

%plot orbit yet another way
figure(3)
hold on
grid on
axis equal
plot3(xvals./R_e, yvals./R_e, zvals./R_e)
[xe, ye, ze] = sphere;
thesurf = surf(xe,ye,ze,'EdgeColor','none');
set(thesurf,'FaceColor',[0 0 1], 'FaceAlpha',0.5,'EdgeAlpha', 0);
title('Molniya Orbit about Earth with J2 Perturbation')
xlabel('x [Earth radii]')
ylabel('y [Earth radii]')
zlabel('z [Earth radii]')
view(-50,15)

%compute orbital elements at each time step
for i = 1:length(trajectory)
    r_curr = trajectory(i,1:3);
    v_curr = trajectory(i,4:6);
    [avals(i), evals(i), ivals(i), RAANvals(i), wvals(i), tavals(i)] = rv_to_oe(r_curr, v_curr, mu_e);
end


%plot orbital elements over time
figure(4)
subplot(6,1,1)
hold on
grid on
box on
plot(tdays,avals./R_e) %semimajor axis (a)
title('Semimajor Axis [Earth radii]')
xlabel('t [days]')

subplot(6,1,2)
hold on
grid on
box on
plot(tdays,evals) %eccentricity (e)
title('Eccentricity')
xlabel('t [days]')

subplot(6,1,3)
hold on
grid on
box on
plot(tdays, ivals) %inclination (i)
title(['Inclination [' char(176) ']'])
xlabel('t [days]')

subplot(6,1,4)
hold on
grid on
box on
plot(tdays, RAANvals) %RAAN (omega)
title(['Right Ascension [' char(176) ']'])
xlabel('t [days]')

subplot(6,1,5)
hold on
grid on
box on
plot(tdays, wvals) %Argument of periapse (w)
title(['Argument of Periapse [' char(176) ']'])
xlabel('t [days]')

subplot(6,1,6)
hold on
grid on
box on
plot(tdays,tavals) %True anomaly (f)
title(['True Anomaly [' char(176) ']'])
xlabel('t [days]')


% % plot orbital elements indiviudally (uncomment if desired)
% %semimajor axis (a)
% figure(5)
% hold on
% grid on
% box on
% plot(tdays, avals)
% title('Semimajor Axis over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('a [km]')
% 
% %eccentricity (e)
% figure(6)
% hold on
% grid on
% box on
% plot(tdays, evals)
% title('Eccentricity over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('e')
% 
% %inclination (i)
% figure(7)
% hold on
% grid on
% box on
% plot(tdays, ivals)
% title('Inclination over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('i [deg]')
% 
% %RAAN (omega)
% figure(8)
% hold on
% grid on
% box on
% plot(tdays, RAANvals)
% title('RAAN over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('omega [deg]')
% 
% %Argument of periapse (w)
% figure(9)
% hold on
% grid on
% box on
% plot(tdays, wvals)
% title('Argument of Periapse over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('w [deg]')
% 
% %True anomaly (f)
% figure(10)
% hold on
% grid on
% box on
% plot(tdays, tavals)
% title('True Anomaly over time with Perturbed Orbit')
% xlabel('t [days]')
% ylabel('f [deg]')
