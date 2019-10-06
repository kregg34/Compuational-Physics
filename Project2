%This program runs the restricted 3 body problem

clear;

%set the total time and step size
totalTime = 10;
steps = 100000;
delta = totalTime/steps;

%constants
M1 = 5.972 * 10^24; %Earth's mass
M2 = 7.35 * 10^22;  %Moon's mass
Mtotal = M1 + M2;
mu2 = M2/Mtotal;
mu2 = 0.03;
mu1 = 1 - mu2;

%dependent varibles
xArr = zeros(steps, 1);
yArr = zeros(steps, 1);
aArr = zeros(steps, 1);
bArr = zeros(steps, 1);

%independent varible
t = zeros(steps, 1);

%initial conditions
initialX = -0.9;
initialY = 0;
initialVx = 0;
initialVy = 0;
initialTime = 0;

%set initial conditions
xArr(1) = initialX;
yArr(1) = initialY;
aArr(1) = initialVx;
bArr(1) = initialVy;
t(1) = initialTime;

%The helper functions
f1 = @(t,x,y,a,b) a; %x dot
f2 = @(t,x,y,a,b) b; %y dot
f3 = @(t,x,y,a,b) x + 2*b - ((1-mu2)*(x+mu2))/(((x+mu2)^2 + y^2)^(3/2)) - (mu2*(x+mu2-1))/(((x-1+mu2)^2 + y^2)^(3/2)); %a dot
f4 = @(t,x,y,a,b) y - 2*a - ((1-mu2)*y)/(((x+mu2)^2 + y^2)^(3/2)) - (mu2*y)/(((x-1+mu2)^2 + y^2)^(3/2)); %b dot

%RK4 algorithm
for n = 1:steps-1
    Kx1 = delta * f1(t(n), xArr(n), yArr(n), aArr(n), bArr(n));
    Ky1 = delta * f2(t(n), xArr(n), yArr(n), aArr(n), bArr(n));
    Ka1 = delta * f3(t(n), xArr(n), yArr(n), aArr(n), bArr(n));
    Kb1 = delta * f4(t(n), xArr(n), yArr(n), aArr(n), bArr(n));

    Kx2 = delta * f1(t(n) + delta/2, xArr(n) + Kx1/2, yArr(n) + Ky1/2, aArr(n) + Ka1/2, bArr(n) + Kb1/2);
    Ky2 = delta * f2(t(n) + delta/2, xArr(n) + Kx1/2, yArr(n) + Ky1/2, aArr(n) + Ka1/2, bArr(n) + Kb1/2);
    Ka2 = delta * f3(t(n) + delta/2, xArr(n) + Kx1/2, yArr(n) + Ky1/2, aArr(n) + Ka1/2, bArr(n) + Kb1/2);
    Kb2 = delta * f4(t(n) + delta/2, xArr(n) + Kx1/2, yArr(n) + Ky1/2, aArr(n) + Ka1/2, bArr(n) + Kb1/2);
    
    Kx3 = delta * f1(t(n) + delta/2, xArr(n) + Kx2/2, yArr(n) + Ky2/2, aArr(n) + Ka2/2, bArr(n) + Kb2/2);
    Ky3 = delta * f2(t(n) + delta/2, xArr(n) + Kx2/2, yArr(n) + Ky2/2, aArr(n) + Ka2/2, bArr(n) + Kb2/2);
    Ka3 = delta * f3(t(n) + delta/2, xArr(n) + Kx2/2, yArr(n) + Ky2/2, aArr(n) + Ka2/2, bArr(n) + Kb2/2);
    Kb3 = delta * f4(t(n) + delta/2, xArr(n) + Kx2/2, yArr(n) + Ky2/2, aArr(n) + Ka2/2, bArr(n) + Kb2/2);
    
    Kx4 = delta * f1(t(n) + delta, xArr(n) + Kx3, yArr(n) + Ky3, aArr(n) + Ka3, bArr(n) + Kb3);
    Ky4 = delta * f2(t(n) + delta, xArr(n) + Kx3, yArr(n) + Ky3, aArr(n) + Ka3, bArr(n) + Kb3);
    Ka4 = delta * f3(t(n) + delta, xArr(n) + Kx3, yArr(n) + Ky3, aArr(n) + Ka3, bArr(n) + Kb3);
    Kb4 = delta * f4(t(n) + delta, xArr(n) + Kx3, yArr(n) + Ky3, aArr(n) + Ka3, bArr(n) + Kb3);

    xArr(n+1) = xArr(n) + (1/6 * (Kx1 + 2*Kx2 + 2*Kx3 + Kx4));
    yArr(n+1) = yArr(n) + (1/6 * (Ky1 + 2*Ky2 + 2*Ky3 + Ky4));
    aArr(n+1) = aArr(n) + (1/6 * (Ka1 + 2*Ka2 + 2*Ka3 + Ka4));
    bArr(n+1) = bArr(n) + (1/6 * (Kb1 + 2*Kb2 + 2*Kb3 + Kb4));
    t(n+1) = t(n) + delta;
end

%constants for zero velocity surface for i'th Lagrange point
C1 = 3 + (3^(4/3))*(mu2^(2/3)) - 10*mu2/3;
C2 = 3 + (3^(4/3))*(mu2^(2/3)) - 14*mu2/3;
C3 = 3 + mu2;
C4 = 3 - mu2;
C5 = 3 - mu2;

%potential function
V1 = @(x,y) 2*mu1/sqrt((x+mu2)^2 + y^2) + 2*mu2/sqrt((x-mu1)^2 + y^2) + x^2 + y^2 - C1;
V2 = @(x,y) 2*mu1/sqrt((x+mu2)^2 + y^2) + 2*mu2/sqrt((x-mu1)^2 + y^2) + x^2 + y^2 - C2;
V3 = @(x,y) 2*mu1/sqrt((x+mu2)^2 + y^2) + 2*mu2/sqrt((x-mu1)^2 + y^2) + x^2 + y^2 - C3;
V4 = @(x,y) 2*mu1/sqrt((x+mu2)^2 + y^2) + 2*mu2/sqrt((x-mu1)^2 + y^2) + x^2 + y^2 - C4;
V5 = @(x,y) 2*mu1/sqrt((x+mu2)^2 + y^2) + 2*mu2/sqrt((x-mu1)^2 + y^2) + x^2 + y^2 - C5;

%Plot everything -------------------------------
figure(1);
xlabel('X')
ylabel('Y')
title('Three Body Problem (Restricted to xy plane)')

%Trajectory
plot(xArr, yArr, 'r');

hold on

%start location
plot(initialX, initialY, 'r*');

%zero velocity curves
fimplicit(V1, 'b');
fimplicit(V2, 'g');
fimplicit(V3, 'm');
fimplicit(V4, 'k');

%L4 and L5 Lagrange points
plot(1/2 - mu2, sqrt(3)/2, 'k*');
plot(1/2 - mu2, -sqrt(3)/2, 'k*');

%planet locations
plot(-mu2, 0, 'bo', 'MarkerFaceColor', 'g');
plot(mu1, 0, 'ko');

hold off
legend('Trajectory', 'Initial position', 'C1', 'C2', 'C3');
