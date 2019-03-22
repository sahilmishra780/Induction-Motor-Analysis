clc
close all

%Induction Motor frameling
global R        %Resistance Matrix for QD framel 
global L        %Inductance Matrix for QD framel
global v        %Given RMS Line-to-Line Voltage
global f        %Frequency in Hz of input 3-phase voltage 
global p        %Number of poles in the Induction Motor 
global J        %Motor Inertia
global Tl       %Load Torque
global k        %Torque Constant
global frame    %Frame of Reference for DQ Model: 1:Synchronous Frame, 2:Stationary frame, 3:Rotor Frame

%Specify Frame of Reference
frame = 3;

%Mechanical Parameters
J = 3.1; Tl = 10; 

%Electrical Parameters
v = 460; f = 60; p = 4; 

%Stator and Rotor Parameters 
rs = 0.0149;                    %Stator resistance
lls = 0.000303;                 %Stator leakage inductance
rr = 0.0093;                    %Rotor resistance
llr = 0.000303;                 %Rotor leakage inductance
lm = 0.0105;                    %Magnetizing inductance
ls = lls + lm;                  %Equivalent stator inductance
lr = llr + lm;                  %Equivalent rotor inductance
sigma = 1 - (lm*lm)/(lr*ls);    %Electromagnetic Torque constant
k = lm/(sigma*ls*lr);           %Multiplier for Torque Constant

%Matrices
R = zeros(4,4); R(1,1) = rs; R(2,2) = rs; R(3,3) = rr; R(4,4) = rr;
L = zeros(4,4); L(1,1) = ls; L(2,2) = ls; L(3,3) = lr; L(4,4) = lr;
L(1,3) = lm; L(2,4) = lm; L(3,1) = lm; L(4,2)= lm;

%4th Order Runge Kutta Solution
tbegin = 0.0;                   %Time of start
tfinal = 4;                     %Time of end
stepsize = 50000;               %Step Size for Runge Kutta
InitialConditions= zeros(6,1);  %Initial Condtions in the order of Lamda_ds, Lamda_qs, Lamda_dr,Lamda_qr, omega_r, theta_r
fh = @Diffsolver;               

%Call the myrk Runge Kutta function
[tspan,x] = myrk(fh,tbegin,tfinal,InitialConditions,stepsize);
x = x';

%Idq Current Vector
idq = zeros(stepsize,4);        

for i=1:1:stepsize
    idq(i,1:4) = L\x(i,1:4)';       %Contains dq currents in order of ids, iqs, idr, iqr
end

%Electromagnetic Torque
etorque = 0.75*p*lm*(idq(:,2).*idq(:,3)-idq(:,4).*idq(:,1));

%Phase current calculation
if frame == 1
    theta = 2*pi*f*tspan';
end
if frame == 2
    theta = zeros(stepsize,1);
end
if frame == 3
    theta = x(:,6);
end

iqds = complex(idq(:,2),-idq(:,1));     %Complex representation of iqds current
iqdr = complex(idq(:,4),-idq(:,3));     %Complex representation of iqdr current
a = exp(1i*(2*pi/3));
b = exp(1i*(-2*pi/3));
theta_r = x(:,6);

%Obtaining Stator and Rotor Phase Currents by Complex Analysis
ias = real(iqds.*exp(1i.*theta));
ibs = real(b*iqds.*exp(1i.*theta));
ics = real(a*iqds.*exp(1i.*theta));
iar = real(iqdr.*exp(1i.*(theta-theta_r)));
ibr = real(b*iqdr.*exp(1i.*(theta-theta_r)));
icr = real(a*iqdr.*exp(1i.*(theta-theta_r)));

%Torque and Speed Plot
subplot(2,1,1)
plot(tspan,etorque)
title('Torque')
xlabel('seconds')
ylabel('N.m')
subplot(2,1,2)
plot(tspan,(60/(pi*p))*x(:,5))
title('Speed')
xlabel('seconds')
ylabel('rpm')

%Stator Phase Currents Plot
figure
subplot(3,1,1)
plot(tspan,ias)
title('ias')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])
subplot(3,1,2)
plot(tspan,ibs)
title('ibs')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])
subplot(3,1,3)
plot(tspan,ics)
title('ics')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])

%Rotor Phase Current Plots
figure
subplot(3,1,1)
plot(tspan,iar)
title('iar')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])
subplot(3,1,2)
plot(tspan,ibr)
title('ibr')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])
subplot(3,1,3)
plot(tspan,icr)
title('icr')
xlabel('seconds')
ylabel('Amps')
axis([tbegin tfinal -2000 2000])
