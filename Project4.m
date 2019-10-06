%Project 4 - Solving the diffusion equation in both 1D and 2D

clearvars;
clc;

%bounds of the varibles
totalT = 1;
totalX = 1;

%boundary conditions
xAtStart = 0;
xAtEnd = 1;

%Step size and number of steps
deltaX = 1/100;
deltaT = (deltaX^2)/2; %value based on stability condition of forward Euler
stepsX = totalX/deltaX;
stepsT = totalT/deltaT;

%Constant, used in all three of the 1D algorithms
global A;
global timeForAnalytic;
global timeStepToStopAt;
timeStepToStopAt = 20;
timeForAnalytic = totalT/timeStepToStopAt;
A = deltaT / (deltaX^2);

%{
initialize the "next position (j+1)" vector
The first and last element are the boundary conditions
%}
nextStepVectorEE = zeros(stepsX+1, 1);
nextStepVectorEE(1) = xAtStart;
nextStepVectorEE(stepsX+1) = xAtEnd;

%initialize the position vector (used for plotting)
xPositionVector = zeros(stepsX + 1, 1);
currentX = 0;
for j = 1:stepsX + 1
    xPositionVector(j) = currentX;
    currentX = currentX + deltaX;
end

%Analytical solution
analyticFunction1D = @(x) x - ((2*sin(pi()*x)*exp(-(pi()^2)*timeForAnalytic))/pi()) + ((2*sin(2*pi()*x)*exp(-4*(pi()^2)*timeForAnalytic))/(2*pi())) - ((2*sin(3*pi()*x)*exp(-9*(pi()^2)*timeForAnalytic))/(3*pi())) + ((2*sin(4*pi()*x)*exp(-16*(pi()^2)*timeForAnalytic))/(4*pi()));

%Explicit forward Euler algorithm-------------------------------------
for t = 1:stepsT
    for n = 2:stepsX
        nextStepVectorEE(n) = A*(nextStepVectorEE(n+1) + nextStepVectorEE(n-1)) + (1-2*A)*(nextStepVectorEE(n));
    end

    %Plot the explicit forward Euler after a few iterations of the method 
    if t == floor(stepsT/timeStepToStopAt)
        figure(1);
        plot(xPositionVector, nextStepVectorEE);
        title('Far from Steady State (with delta X = 1/100)')
        xlabel('X position')
        ylabel('Temperature')
        grid on;
        hold on;
    end
end

%Implicit backward Euler algorithm-----------------------------------

%{
Reset the "next position (j+1)" vector
The first and last element are the boundary conditions
%}
nextStepVectorIE = zeros(stepsX+1, 1);
nextStepVectorIE(1) = xAtStart;
nextStepVectorIE(stepsX+1) = xAtEnd;

%Create the tridiagonal matrix
eulerTridiagonalMatrix = zeros(stepsX-1, stepsX-1);
for i = 1:stepsX-1
    for j = 1:stepsX-1
        if i == j
            eulerTridiagonalMatrix(i,j) = 1+(2*A);
            continue;
        end
 
        if i == j+1 || i == j-1
            eulerTridiagonalMatrix(i,j) = -A;
        end
    end
end

%Part of the calculation for the next time step
addVector = zeros(stepsX-1,1);
addVector(stepsX-1) = A;

%iterate through the time steps
for t = 1:stepsT
    nextStepVectorIE(2:stepsX, 1) = eulerTridiagonalMatrix \ (nextStepVectorIE(2:stepsX, 1) + addVector);
    
    %Plot the implicit backward Euler after a few iterations of the method 
    if t == floor(stepsT/timeStepToStopAt)
        plot(xPositionVector, nextStepVectorIE);
    end
end

%Implicit Crank-Nicolson algorithm-----------------------------------

%{
Reset the "next position (j+1)" vector
The first and last element are the boundary conditions
%}
nextStepVectorCN = zeros(stepsX+1, 1);
nextStepVectorCN(1) = xAtStart;
nextStepVectorCN(stepsX+1) = xAtEnd;

%Create the tridiagonal matrix
cronkTridiagonalMatrix = zeros(stepsX-1, stepsX-1);
for i = 1:stepsX-1
    for j = 1:stepsX-1
        if i == j
            cronkTridiagonalMatrix(i,j) = 1+A;
            continue;
        end
 
        if i == j+1 || i == j-1
            cronkTridiagonalMatrix(i,j) = -A/2;
        end
    end
end

%iterate through the time steps
for t = 1:stepsT
    
    B = (1-A) * nextStepVectorCN(2:stepsX, 1);
    
    C = (A/2) * nextStepVectorCN(3:stepsX, 1);
    C(stepsX-1) = 0;
    
    D(1) = 0;
    D(2:stepsX-1, 1) = (A/2) * nextStepVectorCN(2:stepsX-1, 1);

    nextStepVectorCN(2:stepsX, 1) = cronkTridiagonalMatrix \ (B + C + D + addVector);
    
    %Plot the Crank-Nicolson after a few iterations of the method 
    if t == floor(stepsT/timeStepToStopAt)
        plot(xPositionVector, nextStepVectorCN);
        fplot(@(x) analyticFunction1D(x), [0,1]);
    end
end

legend('I.E','E.E','C.N', 'Analytic');


%2D problem -------------------------------------------

%Total time and space
totalT = 1;
totalX = 1;
totalY = 1;

%Set step sizes and some constants to be used in calculations
deltaX = 1/50;
deltaY = deltaX;
deltaT = (deltaX^2 * deltaY^2)/(2*(deltaX^2 + deltaY^2)); %value based on stability condition of forward Euler

stepsX = totalX/deltaX;
stepsY = totalY/deltaY;
stepsT = totalT/deltaT;

alpha = deltaT/(deltaX^2);
beta = deltaT/(deltaY^2);

%used at each time step in the calculations
nextPositionMatrix = zeros(stepsX+1, stepsY+1);
alphaMatrix = zeros(stepsX-1, stepsY-1);
betaMatrix = zeros(stepsX-1, stepsY-1);

%initialize the position vector (used for plotting)
xPositionVector = zeros(stepsX + 1, 1);
currentX = 0;
for j = 1:stepsX + 1
    xPositionVector(j) = currentX;
    currentX = currentX + deltaX;
end

%initialize the nextPositionMatrix at t=0 (=sin(pi*x)sin(pi*y) at t=0)
for i = 2:stepsX
    for j = 2:stepsY
        nextPositionMatrix(i,j) = sin(pi()*xPositionVector(i)) * sin(pi()*xPositionVector(j));
    end
end

%Plot the 2D matrix
figure(2);
surf(xPositionVector, xPositionVector, nextPositionMatrix)
title('Initial conditions')
axis([0 1 0 1])
zlim([0,1])
view(100,0)

%iterate through the time steps
for t = 1:stepsT+1
    %calculate the alpha and beta matrices
    for i = 2:stepsX
         for j = 2:stepsY
             alphaMatrix(i-1,j-1) = nextPositionMatrix(i+1, j) + nextPositionMatrix(i-1, j);
             betaMatrix(i-1,j-1) = nextPositionMatrix(i, j+1) + nextPositionMatrix(i, j-1);
         end
    end
    
    nextPositionMatrix(2:stepsX, 2:stepsY) = (1-2*alpha-2*beta)*nextPositionMatrix(2:stepsX, 2:stepsY) + alpha*alphaMatrix + beta*betaMatrix;
end

%Plot the 2D matrix
figure(3);
surf(xPositionVector, xPositionVector, nextPositionMatrix)
axis([0 1 0 1])
zlim([0,1])
title('Numeric 2D Solution')
view(100,0)

%Analytical solution
figure(4)
syms x y
analytic2D = sin(pi()*x) * sin(pi()*y) * exp(-2*pi()^2*totalT);
fsurf(analytic2D, [0 1 0 1])
zlim([0,1])
title('Analytic 2D Solution')
view(100,0)
