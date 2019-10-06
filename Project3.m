clear;

%set the total time and step size
numOfPeriods = 300;
totalTime = 2*pi*numOfPeriods;
steps = numOfPeriods * 100;
delta = totalTime/steps;

%global constants
global drivingFrequency;
global drivingAmplitude;
global Q;

%Energy function for a pendulum
energyFunc = @(omega, theta) (1 - cos(theta) + ((omega^2)/2));

%initial conditions and values
length = 1;
gravity = 1;
mass = 1;
dampingConstant = 2;
W_o = sqrt(gravity/length);

initialAngle = 0;
initialAnglarSpeed = 0;
initialTime = 0;

%Constants used in the 2 functions below
drivingFrequency = 2/3;
drivingAmplitude = 1.5;
Q = mass*gravity/(W_o*dampingConstant);

%Dependent varibles
angleArray = zeros(steps, 1);
angularSpeedArray = zeros(steps, 1);
analyticAngle = zeros(steps, 1);
analyticSpeed = zeros(steps, 1);
energyArray = zeros(steps, 1);

%Independent varible
time = zeros(steps, 1);

%Set initial conditions
angleArray(1) = initialAngle;
angularSpeedArray(1) = initialAnglarSpeed;
analyticAngle(1) = initialAngle;
analyticSpeed(1) = initialAnglarSpeed;
time(1) = initialTime;
energyArray(1) = energyFunc(initialAnglarSpeed, initialAngle);

%Functions to be integrated using RK4
f1 = @(theta,thetaDot, t) thetaDot;
f2 = @(theta,thetaDot, t) ((-1)*thetaDot/Q) - sin(theta) + drivingAmplitude*cos(drivingFrequency*t);

%Help functions used below for analytic calculations
minusSq = @() (1 - (drivingFrequency^2));
wOverQ = @() drivingFrequency/Q;
bottom = @() minusSq()^2 + wOverQ()^2;
arg = @(t) drivingFrequency * t;

%Analytic solutions for long t and small angle approx
analyticPos = @(t) drivingAmplitude * ((minusSq()*cos(arg(t)) + wOverQ()*sin(arg(t))) / bottom());
analyticVel = @(t) drivingAmplitude * drivingFrequency * (((-1)*minusSq()*sin(arg(t)) + wOverQ()*cos(arg(t))) / bottom() );

%RK4 algorithm
for n = 1:steps-1

    K1f1 = delta * f1(angleArray(n), angularSpeedArray(n), time(n));
    K1f2 = delta * f2(angleArray(n), angularSpeedArray(n), time(n));

    K2f1 = delta * f1(angleArray(n) + K1f1/2, angularSpeedArray(n) + K1f2/2, time(n) + delta/2);
    K2f2 = delta * f2(angleArray(n) + K1f1/2, angularSpeedArray(n) + K1f2/2, time(n) + delta/2);

    K3f1 = delta * f1(angleArray(n) + K2f1/2, angularSpeedArray(n) + K2f2/2, time(n) + delta/2);
    K3f2 = delta * f2(angleArray(n) + K2f1/2, angularSpeedArray(n) + K2f2/2, time(n) + delta/2);

    K4f1 = delta * f1(angleArray(n) + K3f1, angularSpeedArray(n) + K3f2, time(n) + delta);
    K4f2 = delta * f2(angleArray(n) + K3f1, angularSpeedArray(n) + K3f2, time(n) + delta);

    angleArray(n+1)        =        angleArray(n) + (1/6 * (K1f1 + 2*K2f1 + 2*K3f1 + K4f1));
    angularSpeedArray(n+1) = angularSpeedArray(n) + (1/6 * (K1f2 + 2*K2f2 + 2*K3f2 + K4f2));

    analyticAngle(n+1) = analyticPos(time(n));
    analyticSpeed(n+1) = analyticVel(time(n));

    energyArray(n+1) = energyFunc(angularSpeedArray(n+1), angleArray(n+1));

    %increment the time by one step
    time(n+1) = time(n) + delta;
end

%Create Poincare maps -------------------
nMax = floor(drivingFrequency*totalTime/(2*pi));
poincareTimes = zeros(1, nMax);

for n = 1:nMax
    poincareTimes(n) = (2*pi*n)/drivingFrequency;
end

anglesToPlot = zeros(1, nMax);
speedsToPlot = zeros(1, nMax);
toleranceVal = delta/8;

dim = size(time);
max = dim(1);

for i = 1:max
    for n = 1:nMax
        if (abs(time(i) - poincareTimes(n)) < toleranceVal)
            anglesToPlot(n) = angleArray(i);
            speedsToPlot(n) = angularSpeedArray(i);
        end
    end
end

%Skip the first n values (to avoid transient terms messing with things)
n = 11;
anglesToPlot = anglesToPlot(:, n:nMax);
speedsToPlot = speedsToPlot(:, n:nMax);

%Plot stuff
figure(1);
plot(angularSpeedArray, angleArray);
grid on;
title('Pendulum (A=0.5)')
xlabel('Angular Speed (Radians/time)')
ylabel('Angle (Radians)')

figure(2);
plot(time, angleArray);
grid on;
title('Time vs Angle')
xlabel('Time (unitless)')
ylabel('Angle (Radians)')

figure(3);
scatter(anglesToPlot, speedsToPlot);
grid on;
title('Poincare Map (A=0.5)')
xlabel('Angle (Radians)')
ylabel('Angular Speed (Radians/time)')
