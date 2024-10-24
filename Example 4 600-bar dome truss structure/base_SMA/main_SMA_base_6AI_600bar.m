% Artemisinin slime mold algorithm (ASMA)
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________           

clear all 
close all
clc

disp('The ASMA is tracking the problem');
N=25; % Number of search agents
Function_name='F1'; % Name of the test function, range from F1-F13
MaxIT=500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb, ub, dim, fobj] = Get_Functions_details_Dome600(Function_name);

Times = 2; % Number of independent times you want to run the SMA, LSMA, CLSMA
display(['Number of independent runs: ', num2str(Times)]);
Convergence_curve_ASMA_all = zeros(Times, MaxIT);
Convergence_curve_AO_all = zeros(Times, MaxIT);
Convergence_curve_SMA_all = zeros(Times, MaxIT);
Convergence_curve_LSMA_all = zeros(Times, MaxIT);
Convergence_curve_ESMA_all = zeros(Times, MaxIT);
Convergence_curve_AOSMA_all = zeros(Times, MaxIT);
best_chaos_index = 1; 
chaos_selection_method = 1; % 1 = random per iteration, 2 = best from previous run

disp('The ASMA is tracking the problem');
for i = 1:Times
    %[bestfitness,Leader_pos,Convergence_curve]=ASMA(N,MaxFEs,lb,ub,dim,fobj)
    [Destination_fitness_ASMA(i), bestPositions_ASMA(i,:), Convergence_curve_ASMA(i,:)] =ASMA(N, 12500, lb, ub, dim, fobj);
    Convergence_curve_ASMA_all(i,:) = Convergence_curve_ASMA(i,:);
    display(['The optimal fitness of ASMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_ASMA(i))]);
end

disp('The AO is tracking the problem');
for i = 1:Times
    %[bestfitness,Leader_pos,Convergence_curve]=ASMA(N,MaxFEs,lb,ub,dim,fobj)
    [Destination_fitness_AO(i), bestPositions_AO(i,:), Convergence_curve_AO(i,:)] =AO(N, 12500, lb, ub, dim, fobj);
    Convergence_curve_AO_all(i,:) = Convergence_curve_AO(i,:);
    display(['The optimal fitness of AO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_AO(i))]);
end

disp('The SMA is tracking the problem');
for i = 1:Times
    %[bestfitness,Leader_pos,Convergence_curve]=ASMA(N,MaxFEs,lb,ub,dim,fobj)
    [Destination_fitness_SMA(i), bestPositions_SMA(i,:), Convergence_curve_SMA(i,:)] =SMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_SMA_all(i,:) = Convergence_curve_SMA(i,:);
    display(['The optimal fitness of SMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_SMA(i))]);
end

disp('The LSMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_LSMA(i), bestPositions_LSMA(i,:), Convergence_curve_LSMA(i,:)] = LSMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_LSMA_all(i,:) = Convergence_curve_LSMA(i,:);
    display(['The optimal fitness of LSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_LSMA(i))]);
end

disp('The ESMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_ESMA(i), bestPositions_ESMA(i,:), Convergence_curve_ESMA(i,:)] =ESMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_ESMA_all(i,:) = Convergence_curve_ESMA(i,:);
    display(['The optimal fitness of ESMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_ESMA(i))]);
end


disp('The AOSMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_AOSMA(i), bestPositions_AOSMA(i,:), Convergence_curve_AOSMA(i,:)] = AOSMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_AOSMA_all(i,:) = Convergence_curve_AOSMA(i,:);
    display(['The optimal fitness of AOSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_AOSMA(i))]);
end

[bestfitness_ASMA, index_ASMA] = min(Destination_fitness_ASMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of ASMA is: ', num2str(bestfitness_ASMA)]);
display(['The average fitness of ASMA is: ', num2str(mean(Destination_fitness_ASMA))]);
display(['The standard deviation fitness of ASMA is: ', num2str(std(Destination_fitness_ASMA))]);
display(['The best location of ASMA is: ', num2str(bestPositions_ASMA(index_ASMA,:))]);

[bestfitness_AO, index_AO] = min(Destination_fitness_AO);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of AO is: ', num2str(bestfitness_AO)]);
display(['The average fitness of AO is: ', num2str(mean(Destination_fitness_AO))]);
display(['The standard deviation fitness of AO is: ', num2str(std(Destination_fitness_AO))]);
display(['The best location of AO is: ', num2str(bestPositions_AO(index_AO,:))]);

[bestfitness_SMA, index_SMA] = min(Destination_fitness_SMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of SMA is: ', num2str(bestfitness_SMA)]);
display(['The average fitness of SMA is: ', num2str(mean(Destination_fitness_SMA))]);
display(['The standard deviation fitness of SMA is: ', num2str(std(Destination_fitness_SMA))]);
display(['The best location of SMA is: ', num2str(bestPositions_SMA(index_SMA,:))]);

[bestfitness_LSMA, index_LSMA] = min(Destination_fitness_LSMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of LSMA is: ', num2str(bestfitness_LSMA)]);
display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness_LSMA))]);
display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness_LSMA))]);
display(['The best location of LSMA is: ', num2str(bestPositions_LSMA(index_LSMA,:))]);

[bestfitness_ESMA, index_ESMA] = min(Destination_fitness_ESMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of ESMA is: ', num2str(bestfitness_ESMA)]);
display(['The average fitness of ESMA is: ', num2str(mean(Destination_fitness_ESMA))]);
display(['The standard deviation fitness of ESMA is: ', num2str(std(Destination_fitness_ESMA))]);
display(['The best location of ESMA is: ', num2str(bestPositions_ESMA(index_ESMA,:))]);

[bestfitness_AOSMA, index_AOSMA] = min(Destination_fitness_AOSMA);
disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
display(['The best fitness of AOSMA is: ', num2str(bestfitness_AOSMA)]);
display(['The average fitness of AOSMA is: ', num2str(mean(Destination_fitness_AOSMA))]);
display(['The standard deviation fitness of AOSMA is: ', num2str(std(Destination_fitness_AOSMA))]);
display(['The best location of AOSMA is: ', num2str(bestPositions_LSMA(index_AOSMA,:))]);


% Plotting  convergence curves
figure;
semilogy(Convergence_curve_AO(index_AO,:), 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;
semilogy(mean(Convergence_curve_AO_all), '--', 'LineWidth', 2.5, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;

semilogy(Convergence_curve_SMA(index_SMA,:), 'LineWidth', 3, 'Color', [0.9290, 0.6940, 0.1250]);
hold on;
semilogy(mean(Convergence_curve_SMA_all), '--', 'LineWidth', 2.5, 'Color', [0.9290, 0.6940, 0.1250]);
hold on;

semilogy(Convergence_curve_LSMA(index_LSMA,:), 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;
semilogy(mean(Convergence_curve_LSMA_all), '--', 'LineWidth', 2.5, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;

semilogy(Convergence_curve_ESMA(index_ESMA,:), 'LineWidth', 3, 'Color', [0.4940, 0.1840, 0.5560]);
hold on;
semilogy(mean(Convergence_curve_ESMA_all), '--', 'LineWidth', 2.5, 'Color', [0.4940, 0.1840, 0.5560]);
hold on;

semilogy(Convergence_curve_LSMA(index_AOSMA,:), 'LineWidth', 3, 'Color', [0.4660, 0.6740, 0.1880]);
hold on;
semilogy(mean(Convergence_curve_AOSMA_all), '--', 'LineWidth', 2.5, 'Color', [0.4660, 0.6740, 0.1880]);
hold on;

semilogy(Convergence_curve_ASMA(index_ASMA,:), 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]);
hold on;
semilogy(mean(Convergence_curve_ASMA_all), '--', 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]);
hold on;
xlabel('Iteration');
ylabel('Weight (kg)');


legend('Best, AO', 'Mean, AO','Best, SMA', 'Mean, SMA', 'Best, LSMA', 'Mean, LSMA', 'Best, ESMA', 'Mean, ESMA', 'Best, AOSMA', 'Mean, AOSMA',  'Best, ASMA', 'Mean, ASMA');
title(['Dome 600 bar']);
box on; axis tight; grid on;
current_ylim = ylim;
ylim([current_ylim(1)*0.85, 30000]);  % Set upper y-limit to 3000
% Save the figure
saveas(gcf, 'convergence_curves_Dome_52.png');

% Save best fitness results for all algorithms
best_fitness_results = [bestfitness_ASMA; bestfitness_AO; bestfitness_SMA; bestfitness_LSMA; bestfitness_ESMA; bestfitness_AOSMA]';
csvwrite('best_fitness_results.csv', best_fitness_results);

% Save best positions for all algorithms (as separate files)
csvwrite('best_positions_ASMA.csv', bestPositions_ASMA);
csvwrite('best_positions_AO.csv', bestPositions_AO);
csvwrite('best_positions_SMA.csv', bestPositions_SMA);
csvwrite('best_positions_LSMA.csv', bestPositions_LSMA);
csvwrite('best_positions_ESMA.csv', bestPositions_ESMA);
csvwrite('best_positions_AOSMA.csv', bestPositions_AOSMA);

% Save convergence curves for all algorithms
csvwrite('convergence_curve_ASMA.csv', Convergence_curve_ASMA_all);
csvwrite('convergence_curve_AO.csv', Convergence_curve_AO_all);
csvwrite('convergence_curve_SMA.csv', Convergence_curve_SMA_all);
csvwrite('convergence_curve_LSMA.csv', Convergence_curve_LSMA_all);
csvwrite('convergence_curve_ESMA.csv', Convergence_curve_ESMA_all);
csvwrite('convergence_curve_AOSMA.csv', Convergence_curve_AOSMA_all);

