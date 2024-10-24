clear all 
close all
clc
disp('The ASMA is tracking the problem');
N=25; % Number of search agents

Function_name='F1'; % Name of the test function, range from F1-F13

T=500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details_Dome120(Function_name);

Times=1; %Number of independent times you want to run the ESMA
disp(['Number of independent runs: ', num2str(Times)]);

for i=1:Times
    % Running the SMA
    [Destination_fitness_SMA(i),bestPositions_SMA(i,:),Convergence_curve_SMA(i,:),Exploration_SMA(i,:),Exploitation_SMA(i,:)]= SMA_Div(N,T,lb,ub,dim,fobj);
    disp(['The optimal fitness of SMA is: ', num2str(Destination_fitness_SMA(i))]);

    % Running the ASMA
    [Destination_fitness_ASMA(i),bestPositions_ASMA(i,:),Convergence_curve_ASMA(i,:),Exploration_ASMA(i,:),Exploitation_ASMA(i,:)]= ASMA_Div(N,12500,lb,ub,dim,fobj);
    disp(['The optimal fitness of ASMA is: ', num2str(Destination_fitness_ASMA(i))]);
end

% Get the best results
[bestfitness_SMA,index_SMA]=min(Destination_fitness_SMA);
[bestfitness_ASMA,index_ASMA]=min(Destination_fitness_ASMA);

% Calculate average exploration and exploitation for display in legend
avgExploration_SMA = mean(Exploration_SMA(index_SMA,:));
avgExploitation_SMA = mean(Exploitation_SMA(index_SMA,:));
avgExploration_ASMA = mean(Exploration_ASMA(index_ASMA,:));
avgExploitation_ASMA = mean(Exploitation_ASMA(index_ASMA,:));

% Plot the Exploration and Exploitation percentages over iterations
figure;

% Plotting for SMA
semilogy(Exploration_SMA(index_SMA,:), '-', 'LineWidth', 3, 'Color', [0.9290, 0.6940, 0.1250]); hold on;
semilogy(Exploitation_SMA(index_SMA,:), '--', 'LineWidth', 3, 'Color', [0.9290, 0.6940, 0.1250]); hold on;

% Plotting for ASMA
semilogy(Exploration_ASMA(index_ASMA,:), '-', 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]); hold on;
semilogy(Exploitation_ASMA(index_ASMA,:), '--', 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]); hold on;

% Set axis limits and turn off grid
ylim([0 100]);
box on; axis tight; grid off;

% Set axis labels and customized legend
xlabel('Iteration');
ylabel('Percentage (%)');
legend(sprintf('Avg Exploration SMA %%=%.2f', avgExploration_SMA), ...
       sprintf('Avg Exploitation SMA %%=%.2f', avgExploitation_SMA), ...
       sprintf('Avg Exploration ASMA %%=%.2f', avgExploration_ASMA), ...
       sprintf('Avg Exploitation ASMA %%=%.2f', avgExploitation_ASMA));
