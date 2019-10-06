%Project 5 - The 2D Ising Model (assuming a square lattice and h=0)
clc;
clear;

%% Defining varibles

startingTemp = 0;
maxTemp = 10;
tempIncrement = 0.05;
latticeSize = 2;
numberOfSpins = latticeSize^2;
numberOfCycles = 1000;
k = 1;
J = 1;

magnetizations = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);
energies = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);
susceptibilities = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);
heatCapacities = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);
acceptedConfigs = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);

temperatures = startingTemp:tempIncrement:ceil(maxTemp);

%% Set up a random spin configuration

spinConfigurationInitial = rand(latticeSize, latticeSize);

for i = 1:latticeSize
    for j = 1:latticeSize
        if spinConfigurationInitial(i,j) > 0.5
            spinConfigurationInitial(i,j) = 1;
        else
            spinConfigurationInitial(i,j) = -1;
        end
    end
end

%% Iterate over temperatures and run the Metropolis algorithm

indexForNewTemps = 1;

for temp = temperatures
    
    %Reset the configuration to the initial one for each new temperature
    spinConfiguration = spinConfigurationInitial;
    
    %Thermodynamic quantities
    expectedE = 0;
    expectedESquared = 0;
    expectedM = 0;
    expectedMSquared = 0;
    accepted = 0;

    %Loop through all the Monte Carlo cycles
    for i = 1:numberOfCycles
        
        %Randomly choice the spins to visit (and visit them all without
        %repeating)
        orderOfSelections = randperm(numberOfSpins, numberOfSpins);
        
        %This for loop is one Monte Carlo cycle (goes through all spins
        %randomly)
        for randomIndex = 1:numberOfSpins
            choice = orderOfSelections(randomIndex);

            col = mod(choice, latticeSize);

            if col == 0
                col = latticeSize;
            end

            row = floor(choice/(latticeSize+0.1)) + 1;
            
            %Get the values of the indexes around the random XY point
            rowPlusOne = mod(row, latticeSize) + 1;
            rowMinusOne  = mod(row - 2, latticeSize) + 1;
            colPlusOne = mod(col, latticeSize) + 1;
            colMinusOne = mod(col - 2, latticeSize) + 1;

            rowPlusOnePoint = spinConfiguration(rowPlusOne, col);
            rowMinusOnePoint = spinConfiguration(rowMinusOne, col);
            colPlusOnePoint = spinConfiguration(row, colPlusOne);
            colMinusOnePoint = spinConfiguration(row, colMinusOne);

            sumOfPoints = rowPlusOnePoint + rowMinusOnePoint + colPlusOnePoint + colMinusOnePoint;
            
            %Energy change due to flipping the spin
            deltaE = 2*J*sumOfPoints*spinConfiguration(row, col);
            
            %Probablility used to determine if we keep the new configuration
            p_acceptance = exp(-deltaE/(k*temp));

            %Decide if we will keep the new spin configuration or not
            if deltaE <= 0 || rand() < p_acceptance
                %New configuation is accepted
                accepted = accepted + 1;
                spinConfiguration(row, col) = -1 * spinConfiguration(row, col);
            end
        end
        
        energy = findTotalEnergy(spinConfiguration, J);
        magnetization = sum(spinConfiguration(:));
        
        expectedE = expectedE + energy;
        expectedESquared = expectedESquared + (energy^2);
        expectedM = expectedM + abs(magnetization);
        expectedMSquared = expectedMSquared + (magnetization^2);
    end
    
    divisor1 = numberOfCycles * numberOfSpins;
    divisor2 = ((numberOfCycles)^2) * numberOfSpins;
    
    energies(indexForNewTemps, 1) = expectedE / divisor1;
    magnetizations(indexForNewTemps, 1) = expectedM / divisor1;
    susceptibilities(indexForNewTemps, 1) = ((expectedMSquared/divisor1) - ((expectedM^2)/divisor2)) / (temp);
    heatCapacities(indexForNewTemps, 1) = ((expectedESquared/divisor1) - ((expectedE^2)/divisor2)) / (temp^2);
    acceptedConfigs(indexForNewTemps, 1) = accepted / divisor1;
    
    indexForNewTemps = indexForNewTemps + 1;
end

%% Analytic Solution for 2x2 case
numOfSpins = 4;

energiesAnalyticTwobyTwo = zeros(ceil((maxTemp - startingTemp)/tempIncrement + 1), 1);
magTwobyTwo = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1);
capTwobyTwo = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1); 
susTwobyTwo = zeros((maxTemp - startingTemp)/tempIncrement + 1, 1); 

twoByTwoEnergy = @(T) (8*(exp(-8/T) - exp(8/T)))  /   (( exp(8/T) +  exp(-8/T) + 6)) / numOfSpins;
twoByTwoMag = @(T) (8*exp(8/T) + 16)  /  ((2*exp(8/T) + 2*exp(-8/T) + 12)) / numOfSpins;
twoByTwoCap = @(T) ((128*exp(8/T))/(T^2)) * ( (2*exp(8/T) + 3*exp(16/T) + 3) / (((6*exp(8/T) + exp(16/T) + 1)^2)) ) / numOfSpins;
twoByTwoSus = @(T) (1/T) * ( ((32*exp(8/T) + 32)/(2*exp(8/T) + 12 + 2*exp(-8/T))) - (((8*exp(8/T) + 16) / (2*exp(8/T) + 12 + 2*exp(-8/T)))^2) ) / numOfSpins;

index = 1;
for temp = temperatures
    energiesAnalyticTwobyTwo(index, 1) = twoByTwoEnergy(temp);
    magTwobyTwo(index, 1) = twoByTwoMag(temp);
    capTwobyTwo(index, 1) = twoByTwoCap(temp);
    susTwobyTwo(index, 1) = twoByTwoSus(temp);
    index = index + 1;
end

%% Plot everything against temperature

figure (1)
plot(temperatures, magnetizations)
hold on;
plot(temperatures, magTwobyTwo);
hold off;
title('Temperature Vs. Magnetization')
xlabel('Temperature')
ylabel('Magnetization')
ylim([-0.1 1.1])
legend('Monte Carlo', 'Analytic')

figure (2)
plot(temperatures, energies)
hold on;
plot(temperatures, energiesAnalyticTwobyTwo, 'r')
hold off;
title('Temperature Vs. Energy')
xlabel('temperature')
ylabel('Energy')
legend('Monte Carlo', 'Analytic', 'Location', 'Southeast')

figure (3)
plot(temperatures, susceptibilities)
hold on;
plot(temperatures, susTwobyTwo)
hold off;
title('Temperature Vs. Susceptibility')
xlabel('Temperature')
ylabel('Susceptibility')
legend('Monte Carlo', 'Analytic', 'Location', 'Southeast')

figure (4)
plot(temperatures, heatCapacities)
hold on;
plot(temperatures, capTwobyTwo, 'r')
hold off;
title('Temperature Vs. Heat Capacity')
xlabel('Temperature')
ylabel('Heat Capacity')
legend('Monte Carlo', 'Analytic')

figure (100)
plot(temperatures, acceptedConfigs)
title('Temperature Vs. Accepted # of Configs')
xlabel('Temperature')
ylabel('Accepted # of Configs')


%% Looking at the 20x20 lattice as a function of Monte Carlo sweeps
clear;

%Some constants
temp = 2.4;
latticeSize = 20;
numberOfSpins = latticeSize^2;

%Inital 20x20 configuations
orderedConfigurationInitial = ones(latticeSize,latticeSize);
randomConfigurationInitial = rand(latticeSize, latticeSize);

for i = 1:latticeSize
    for j = 1:latticeSize
        if randomConfigurationInitial(i,j) > 0.5
            randomConfigurationInitial(i,j) = 1;
        else
            randomConfigurationInitial(i,j) = -1;
        end
    end
end

%Change this line if you want to switch to ordered/unordered
latticeToBeUsedInitial = orderedConfigurationInitial;

J = 1;
k = 1;
minCycle = 5;
maxCycle = 200;
increCycle = 100;
cycles = minCycle:increCycle:maxCycle;

magnetizations = zeros(floor((maxCycle - minCycle)/increCycle + 1), 1);
energies = zeros(floor((maxCycle - minCycle)/increCycle + 1), 1);
acceptedConfigs = zeros(floor((maxCycle - minCycle)/increCycle + 1), 1);

%Go through the Monte Carlo algorithm
indexForNumCycles = 1;

for numberOfCycles = cycles
    
    %Reset the configuration to the initial one for each new temperature
    spinConfiguration = latticeToBeUsedInitial;
    
    %Thermodynamic quantities
    expectedE = 0;
    expectedESquared = 0;
    expectedM = 0;
    expectedMSquared = 0;
    accepted = 0;

    %Loop through all the Monte Carlo cycles
    for i = 1:numberOfCycles
        
        %Randomly choice the spins to visit (and visit them all without
        %repeating)
        orderOfSelections = randperm(numberOfSpins, numberOfSpins);
        
        %This for loop is one Monte Carlo cycle (goes through all spins
        %randomly)
        for randomIndex = 1:numberOfSpins
            choice = orderOfSelections(randomIndex);

            col = mod(choice, latticeSize);

            if col == 0
                col = latticeSize;
            end

            row = floor(choice/(latticeSize+0.1)) + 1;
            
            %Get the values of the indexes around the random XY point
            rowPlusOne = mod(row, latticeSize) + 1;
            rowMinusOne  = mod(row - 2, latticeSize) + 1;
            colPlusOne = mod(col, latticeSize) + 1;
            colMinusOne = mod(col - 2, latticeSize) + 1;

            rowPlusOnePoint = spinConfiguration(rowPlusOne, col);
            rowMinusOnePoint = spinConfiguration(rowMinusOne, col);
            colPlusOnePoint = spinConfiguration(row, colPlusOne);
            colMinusOnePoint = spinConfiguration(row, colMinusOne);

            sumOfPoints = rowPlusOnePoint + rowMinusOnePoint + colPlusOnePoint + colMinusOnePoint;
            
            %Energy change due to flipping the spin
            deltaE = 2*J*sumOfPoints*spinConfiguration(row, col);
            
            %Probablility used to determine if we keep the new configuration
            p_acceptance = exp(-deltaE/(k*temp));

            %Decide if we will keep the new spin configuration or not
            if deltaE <= 0 || rand() < p_acceptance
                %New configuation is accepted
                accepted = accepted + 1;
                spinConfiguration(row, col) = -1 * spinConfiguration(row, col);
            end
        end
        
        energy = findTotalEnergy(spinConfiguration, J);
        magnetization = sum(spinConfiguration(:));
        
        expectedE = expectedE + energy;
        expectedESquared = expectedESquared + (energy^2);
        expectedM = expectedM + abs(magnetization);
        expectedMSquared = expectedMSquared + (magnetization^2);
    end
    
    divisor1 = numberOfCycles * numberOfSpins;
    
    energies(indexForNumCycles, 1) = expectedE / divisor1;
    magnetizations(indexForNumCycles, 1) = expectedM / divisor1;
    acceptedConfigs(indexForNumCycles, 1) = accepted / divisor1;
    
    indexForNumCycles = indexForNumCycles + 1;
end

%% Plot the 20x20 case

tempStr = num2str(temp);

figure (5)
plot(cycles, magnetizations)
title(['Number of cycles Vs. Magnetization, T=', tempStr])
xlabel('Number of Monte Carlo cycles')
ylabel('Magnetization')
ylim([-0.1 1.1])

figure (6)
plot(cycles, energies)
title(['Number of cycles Vs. Energy, T=', tempStr])
xlabel('Number of Monte Carlo cycles')
ylabel('Energy')

figure (7)
plot(cycles, acceptedConfigs)
title(['Number of cycles Vs. Accepted number of configs, T=', tempStr])
xlabel('Number of Monte Carlo cycles')
ylabel('Accepted number of configs')

%% Find the probability density of the energies P(E)

clear;

%Some constants
temp = 2.4;
latticeSize = 20;
numberOfSpins = latticeSize^2;
J = 1;
k = 1;

%Inital 20x20 configuations
spinConfiguration = rand(latticeSize, latticeSize);

for i = 1:latticeSize
    for j = 1:latticeSize
        if spinConfiguration(i,j) > 0.5
            spinConfiguration(i,j) = 1;
        else
            spinConfiguration(i,j) = -1;
        end
    end
end

cyclesTillEq = 10000;

%Bring the spin configuation to an equilibrium point (using 10,000 cycles)
for i = 1:cyclesTillEq

    %Randomly choice the spins to visit (and visit them all without
    %repeating)
    orderOfSelections = randperm(numberOfSpins, numberOfSpins);

    %This for loop is one Monte Carlo cycle (goes through all spins
    %randomly)
    for randomIndex = 1:numberOfSpins
        choice = orderOfSelections(randomIndex);

        col = mod(choice, latticeSize);

        if col == 0
            col = latticeSize;
        end

        row = floor(choice/(latticeSize+0.1)) + 1;

        %Get the values of the indexes around the random XY point
        rowPlusOne = mod(row, latticeSize) + 1;
        rowMinusOne  = mod(row - 2, latticeSize) + 1;
        colPlusOne = mod(col, latticeSize) + 1;
        colMinusOne = mod(col - 2, latticeSize) + 1;

        rowPlusOnePoint = spinConfiguration(rowPlusOne, col);
        rowMinusOnePoint = spinConfiguration(rowMinusOne, col);
        colPlusOnePoint = spinConfiguration(row, colPlusOne);
        colMinusOnePoint = spinConfiguration(row, colMinusOne);

        sumOfPoints = rowPlusOnePoint + rowMinusOnePoint + colPlusOnePoint + colMinusOnePoint;

        %Energy change due to flipping the spin
        deltaE = 2*J*sumOfPoints*spinConfiguration(row, col);

        %Probablility used to determine if we keep the new configuration
        p_acceptance = exp(-deltaE/(k*temp));

        %Decide if we will keep the new spin configuration or not
        if deltaE <= 0 || rand() < p_acceptance
            spinConfiguration(row, col) = -1 * spinConfiguration(row, col);
        end
    end
end

numberOfCycles = 50000;
energies = zeros(numberOfCycles, 1);
index = 1;

%Thermodynamic quantities
expectedE = 0;
expectedESquared = 0;

%Loop through all the Monte Carlo cycles now that the system is in equilibrium
for i = 1:numberOfCycles

    %Randomly choice the spins to visit (and visit them all without
    %repeating)
    orderOfSelections = randperm(numberOfSpins, numberOfSpins);

    %This for loop is one Monte Carlo cycle (goes through all spins
    %randomly)
    for randomIndex = 1:numberOfSpins
        choice = orderOfSelections(randomIndex);

        col = mod(choice, latticeSize);

        if col == 0
            col = latticeSize;
        end

        row = floor(choice/(latticeSize+0.1)) + 1;

        %Get the values of the indexes around the random XY point
        rowPlusOne = mod(row, latticeSize) + 1;
        rowMinusOne  = mod(row - 2, latticeSize) + 1;
        colPlusOne = mod(col, latticeSize) + 1;
        colMinusOne = mod(col - 2, latticeSize) + 1;

        rowPlusOnePoint = spinConfiguration(rowPlusOne, col);
        rowMinusOnePoint = spinConfiguration(rowMinusOne, col);
        colPlusOnePoint = spinConfiguration(row, colPlusOne);
        colMinusOnePoint = spinConfiguration(row, colMinusOne);

        sumOfPoints = rowPlusOnePoint + rowMinusOnePoint + colPlusOnePoint + colMinusOnePoint;

        %Energy change due to flipping the spin
        deltaE = 2*J*sumOfPoints*spinConfiguration(row, col);

        %Probablility used to determine if we keep the new configuration
        p_acceptance = exp(-deltaE/(k*temp));

        %Decide if we will keep the new spin configuration or not
        if deltaE <= 0 || rand() < p_acceptance
            %New configuation is accepted
            spinConfiguration(row, col) = -1 * spinConfiguration(row, col);
        end
    end

    energy = findTotalEnergy(spinConfiguration, J);
    magnetization = sum(spinConfiguration(:));

    energies(index, 1) = energy;
    
    expectedE = expectedE + energy;
    expectedESquared = expectedESquared + (energy^2);
    
    index = index + 1;
end

divisor1 = numberOfCycles * numberOfSpins;
divisor2 = ((numberOfCycles)^2) * numberOfSpins;

varE = ((expectedESquared/divisor1) - ((expectedE^2)/divisor2))
probE = tabulate(energies);
tempStr = num2str(temp);

%% Plot P(E)

figure (8)
plot(probE(:,1), probE(:,3)/100)
title(['Energy Vs. P(E), T=', tempStr])
xlabel('Energy')
ylabel('Probability of a given energy P(E)')
ylim([0 1])

     
%% Calculate the energy of a given configuation

function energy = findTotalEnergy(config, J)

    energyAmount = 0;
    latticeSizeX = size(config,1);
    
    for x = 1:latticeSizeX
        for y = 1:latticeSizeX
            XPlusOne = mod(x, latticeSizeX) + 1;
            XMinusOne  = mod(x - 2, latticeSizeX) + 1;
            YMinusOne = mod(y - 2, latticeSizeX) + 1;
            YPlusOne = mod(y, latticeSizeX) + 1;

            XPlusOnePoint = config(y, XPlusOne);
            XMinusOnePoint = config(y, XMinusOne);
            YPlusOnePoint = config(YPlusOne, x);
            YMinusOnePoint = config(YMinusOne, x);
            
            sumOfPoints = YPlusOnePoint + YMinusOnePoint + XPlusOnePoint + XMinusOnePoint;
            energyAmount = energyAmount + (sumOfPoints*config(y, x));
        end
    end
    
    energy  = -J * energyAmount / 2;
end
