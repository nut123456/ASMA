%Artemisinin slime mold algorithm (ASMA)
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________ 
clearvars 
close all
clc

disp('The algorithm is tracking the problem');
N = 25; % Number of search agents
Function_name = 'F1'; % Name of the test function, range from F1-F13
MaxIT = 500; % Maximum number of iterations

% Load details of the selected benchmark function
[lb, ub, dim, fobj] = Get_Functions_details_F1_F23(Function_name);

Times = 1; % Number of independent times you want to run all algorithms
%display(['Number of independent runs: ', num2str(Times)]);

% Initialize storage for convergence curves
Convergence_curve_SMA_all = zeros(Times, MaxIT);
Convergence_curve_LSMA_all = zeros(Times, MaxIT);
Convergence_curve_CLSMA_all = zeros(Times, MaxIT);

Convergence_curve_ASMA_all = zeros(Times, MaxIT);                          %1
% Add variables for new algorithms (AO, ECO, INFO, etc.)
Convergence_curve_AO_all = zeros(Times, MaxIT);                          %2
Convergence_curve_SMA_all = zeros(Times, MaxIT);                          %3
Convergence_curve_QIO_all = zeros(Times, MaxIT);                          %4
Convergence_curve_PLO_all = zeros(Times, MaxIT);                          %5
Convergence_curve_PO_all = zeros(Times, MaxIT);                          %6
Convergence_curve_ETO_all = zeros(Times, MaxIT);                          %7
% Add variables for base-SMA algorithms
Convergence_curve_ESMA_all = zeros(Times, MaxIT);                          %8
Convergence_curve_AOSMA_all = zeros(Times, MaxIT);                          %9
Convergence_curve_LSMA_all = zeros(Times, MaxIT);                          %10
Convergence_curve_RCLSMAOA_all = zeros(Times, MaxIT);                          %11

best_chaos_index = 1; 
chaos_selection_method = 1; % 1 = random per iteration, 2 = best from previous run

%1🌿🦟 ASMA Optimization Algorithm (2024)

for i = 1:Times
    [Destination_fitness_ASMA, bestPositions_ASMA, Convergence_curve_ASMA] = ASMA(N, 12500 ,lb, ub, dim,fobj ); %ASMA(N,MaxFEs,lb,ub,dim,fobj) 
    Convergence_curve_ASMA_all = Convergence_curve_ASMA;
    %display(['The optimal fitness of AO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_AO)]);
    display(['The best fitness of ASMA is: ', num2str(min(Convergence_curve_ASMA_all))]);
    disp('The  ASMA finish');
end

%2📜 AO 
for i = 1:Times
    [Destination_fitness_AO, bestPositions_AO, Convergence_curve_AO] = AO(fobj, lb, ub, dim,N, 12500 );  %AO(fobj,lb,ub,dim,N,MaxFEs)
    Convergence_curve_AO_all = Convergence_curve_AO;
    %display(['The optimal fitness of ECO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_ECO)]);
    display(['The best fitness of AO is: ', num2str(min(Convergence_curve_AO_all))]);
    disp('The  AO finish');
end

%3 SMA Algorithm
for i = 1:Times
    [Destination_fitness_SMA(i), bestPositions_SMA(i,:), Convergence_curve_SMA(i,:)] = SMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_SMA_all(i,:) = Convergence_curve_SMA(i,:);
    %display(['The optimal fitness of SMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_SMA(i))]);
    display(['The best fitness of SMA is: ', num2str(min(Convergence_curve_SMA_all))]);
    disp('The  SMA finish');
end

% %4 QIO Algorithm (2023)
% for i = 1:Times
%     [Destination_fitness_QIO, bestPositions_QIO, Convergence_curve_QIO] = QIO(N, MaxIT, lb, ub, dim ,fobj); 
%     Convergence_curve_QIO_all = Convergence_curve_QIO;
%     %display(['The optimal fitness of QIO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_QIO(i))]);
%     display(['The best fitness of AO is: ', num2str(Destination_fitness_QIO)]);
%     disp('The  QIO finish');
% end

%5 📜 Polar Lights Optimizer (PLO)
for i = 1:Times
    [ bestPositions_PLO,Destination_fitness_PLO, Convergence_curve_PLO] = PLO(N, 12500, lb, ub, dim,fobj ); %PLO(N,MaxFEs,lb,ub,dim,fobj)
    Convergence_curve_PLO_all = Convergence_curve_PLO;
    %display(['The optimal fitness of PLO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_PLO(i))]);
    display(['The best fitness of PLO is: ', num2str(min(Convergence_curve_PLO_all))]);
    disp('The  PLO finish');
end

%6 Parrot Optimizer (PO) Algorithm
for i = 1:Times
    [Destination_fitness_PO, bestPositions_PO, Convergence_curve_PO] = PO(N, MaxIT, lb, ub, dim ,fobj ); %PO(N, Max_iter, lb, ub, dim, fobj)
    Convergence_curve_PO_all = Convergence_curve_PO;
    %display(['The optimal fitness of PO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_PO(i))]);
    display(['The best fitness of PO is: ', num2str(min(Convergence_curve_PO_all))]);
    disp('The  PO finish');
end

%7 📜 ETO Optimization (2024)
for i = 1:Times
    [Destination_fitness_ETO, bestPositions_ETO, Convergence_curve_ETO] = ETO(N, 500, lb, ub, dim,fobj ); %ETO(N,Max_Iter,LB,UB,Dim,Fobj)
    Convergence_curve_ETO_all = Convergence_curve_ETO;
    %display(['The optimal fitness of ETO (Run ', num2str(i), ') is: ', num2str(Destination_fitness_ETO(i))]);
    display(['The best fitness of ETO is: ', num2str(min(Convergence_curve_ETO_all))]);
    disp('The  ETO finish');
end

%8 ESMA Algorithm
disp('The ESMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_ESMA(i), bestPositions_ESMA(i,:), Convergence_curve_ESMA(i,:)] = ESMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_ESMA_all(i,:) = Convergence_curve_ESMA(i,:);
    %display(['The optimal fitness of ESMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_ESMA(i))]);
    display(['The best fitness of ESMA is: ', num2str(min(Convergence_curve_ESMA_all))]);
    disp('The  ESMA finish');
end

%9 AOSMA Algorithm
disp('The AOSMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_AOSMA(i), bestPositions_AOSMA(i,:), Convergence_curve_AOSMA(i,:)] = AOSMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_AOSMA_all(i,:) = Convergence_curve_AOSMA(i,:);
    %display(['The optimal fitness of AOSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_AOSMA(i))]);
    display(['The best fitness of AOSMA is: ', num2str(min(Convergence_curve_AOSMA_all))]);
    disp('The  AOSMA finish');
end

%10 LSMA Algorithm
disp('The LSMA is tracking the problem');
for i = 1:Times
    [Destination_fitness_LSMA(i), bestPositions_LSMA(i,:), Convergence_curve_LSMA(i,:)] = LSMA(N, MaxIT, lb, ub, dim, fobj);
    Convergence_curve_LSMA_all(i,:) = Convergence_curve_LSMA(i,:);
    %display(['The optimal fitness of LSMA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_LSMA(i))]);
    display(['The best fitness of LSMA is: ', num2str(Destination_fitness_LSMA)]);
    disp('The  LSMA finish');
end

% %11 RCLSMAOA Algorithm
% disp('The RCLSMAOA is tracking the problem');
% for i = 1:Times
%     [Destination_fitness_RCLSMAOA(i), bestPositions_RCLSMAOA(i,:), Convergence_curve_RCLSMAOA(i,:)] = RCLSMAOA(N, MaxIT, lb, ub, dim, fobj);
%     Convergence_curve_RCLSMAOA_all(i,:) = Convergence_curve_RCLSMAOA(i,:);
%     %display(['The optimal fitness of RCLSMAOA (Run ', num2str(i), ') is: ', num2str(Destination_fitness_RCLSMAOA(i))]);
%     display(['The best fitness of RCLSMAOA is: ', num2str(Destination_fitness_RCLSMAOA)]);
%     disp('The  RCLSMAOA finished');
% end





% % --------Best Fitness, Average Fitness, Standard Deviation and Best Solution for each algorithm--------
% 
% [bestfitness_AO, index_AO] = min(Destination_fitness_AO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of AO is: ', num2str(bestfitness_AO)]);
% display(['The average fitness of AO is: ', num2str(mean(Destination_fitness_AO))]);
% display(['The standard deviation fitness of AO is: ', num2str(std(Destination_fitness_AO))]);
% display(['The best location of AO is: ', num2str(bestPositions_AO)]);
% 
% [bestfitness_ECO, index_ECO] = min(Destination_fitness_ECO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of ECO is: ', num2str(bestfitness_ECO)]);
% display(['The average fitness of ECO is: ', num2str(mean(Destination_fitness_ECO))]);
% display(['The standard deviation fitness of ECO is: ', num2str(std(Destination_fitness_ECO))]);
% display(['The best location of ECO is: ', num2str(bestPositions_ECO)]);
% 
% [bestfitness_FATA, index_FATA] = min(Destination_fitness_FATA);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of FATA is: ', num2str(bestfitness_FATA)]);
% display(['The average fitness of FATA is: ', num2str(mean(Destination_fitness_FATA))]);
% display(['The standard deviation fitness of FATA is: ', num2str(std(Destination_fitness_FATA))]);
% display(['The best location of FATA is: ', num2str(bestPositions_FATA)]);
% 
% % [bestfitness_HGS, index_HGS] = min(Destination_fitness_HGS);
% % disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% % display(['The best fitness of HGS is: ', num2str(bestfitness_HGS)]);
% % display(['The average fitness of HGS is: ', num2str(mean(Destination_fitness_HGS))]);
% % display(['The standard deviation fitness of HGS is: ', num2str(std(Destination_fitness_HGS))]);
% % display(['The best location of HGS is: ', num2str(bestPositions_HGS(index_HGS,:))]);
% 
% [bestfitness_MGO, index_MGO] = min(Destination_fitness_MGO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of MGO is: ', num2str(bestfitness_MGO)]);
% display(['The average fitness of MGO is: ', num2str(mean(Destination_fitness_MGO))]);
% display(['The standard deviation fitness of MGO is: ', num2str(std(Destination_fitness_MGO))]);
% display(['The best location of MGO is: ', num2str(bestPositions_MGO)]);
% 
% [bestfitness_PO, index_PO] = min(Destination_fitness_PO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of PO is: ', num2str(bestfitness_PO)]);
% display(['The average fitness of PO is: ', num2str(mean(Destination_fitness_PO))]);
% display(['The standard deviation fitness of PO is: ', num2str(std(Destination_fitness_PO))]);
% display(['The best location of PO is: ', num2str(bestPositions_PO)]);
% 
% [bestfitness_PLO, index_PLO] = min(Destination_fitness_PLO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of PLO is: ', num2str(bestfitness_PLO)]);
% display(['The average fitness of PLO is: ', num2str(mean(Destination_fitness_PLO))]);
% display(['The standard deviation fitness of PLO is: ', num2str(std(Destination_fitness_PLO))]);
% display(['The best location of PLO is: ', num2str(bestPositions_PLO)]);
% 
% [bestfitness_RIME, index_RIME] = min(Destination_fitness_RIME);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of RIME is: ', num2str(bestfitness_RIME)]);
% display(['The average fitness of RIME is: ', num2str(mean(Destination_fitness_RIME))]);
% display(['The standard deviation fitness of RIME is: ', num2str(std(Destination_fitness_RIME))]);
% display(['The best location of RIME is: ', num2str(bestPositions_RIME(index_RIME,:))]);
% 
% [bestfitness_INFO, index_INFO] = min(Destination_fitness_INFO);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of INFO is: ', num2str(bestfitness_INFO)]);
% display(['The average fitness of INFO is: ', num2str(mean(Destination_fitness_INFO))]);
% display(['The standard deviation fitness of INFO is: ', num2str(std(Destination_fitness_INFO))]);
% display(['The best location of INFO is: ', num2str(bestPositions_INFO(index_INFO,:))]);
% % 
% % [bestfitness_HHO_JOS, index_HHO_JOS] = min(Destination_fitness_HHO_JOS);
% % disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% % display(['The best fitness of HHO_JOS is: ', num2str(bestfitness_HHO_JOS)]);
% % display(['The average fitness of HHO_JOS is: ', num2str(mean(Destination_fitness_HHO_JOS))]);
% % display(['The standard deviation fitness of HHO_JOS is: ', num2str(std(Destination_fitness_HHO_JOS))]);
% % display(['The best location of HHO_JOS is: ', num2str(bestPositions_HHO_JOS(index_HHO_JOS,:))]);
% 
% [bestfitness_SMA, index_SMA] = min(Destination_fitness_SMA);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of SMA is: ', num2str(bestfitness_SMA)]);
% display(['The average fitness of SMA is: ', num2str(mean(Destination_fitness_SMA))]);
% display(['The standard deviation fitness of SMA is: ', num2str(std(Destination_fitness_SMA))]);
% display(['The best location of SMA is: ', num2str(bestPositions_SMA(index_SMA,:))]);
% 
% [bestfitness_LSMA, index_LSMA] = min(Destination_fitness_LSMA);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of LSMA is: ', num2str(bestfitness_LSMA)]);
% display(['The average fitness of LSMA is: ', num2str(mean(Destination_fitness_LSMA))]);
% display(['The standard deviation fitness of LSMA is: ', num2str(std(Destination_fitness_LSMA))]);
% display(['The best location of LSMA is: ', num2str(bestPositions_LSMA(index_LSMA,:))]);
% 
% [bestfitness_CLSMA, index_CLSMA] = min(Destination_fitness_CLSMA);
% disp('--------Best Fitness, Average Fitness, Standard Deviation and Best Solution--------');
% display(['The best fitness of CELSMA is: ', num2str(bestfitness_CLSMA)]);
% display(['The average fitness of CELSMA is: ', num2str(mean(Destination_fitness_CLSMA))]);
% display(['The standard deviation fitness of CELSMA is: ', num2str(std(Destination_fitness_CLSMA))]);
% display(['The best location of CELSMA is: ', num2str(bestPositions_LSMA(index_CLSMA,:))]);



% Plotting  convergence curves
figure;
% figure('Position', [269, 240, 460, 290])
% Plot for ASMA
semilogy(Convergence_curve_ASMA_all, '-', 'LineWidth', 2.5, 'Color', [0, 0.4470, 0.7410]); % Blue solid with square markers
hold on;

% Plot for AO
semilogy(Convergence_curve_AO_all, ':', 'LineWidth', 2.5, 'Color', [0.8500, 0.3250, 0.0980]); % Red dotted with diamond markers
hold on;

% Plot for SMA
semilogy(Convergence_curve_SMA_all, '--', 'LineWidth', 2.5, 'Color', [0.9290, 0.6940, 0.1250]); % Yellow dashed with circle markers
hold on;

% Plot for QIO
semilogy(Convergence_curve_QIO_all, '-', 'LineWidth', 2.5, 'Color', [0.4940, 0.1840, 0.5560]); % Purple solid with triangle-up markers
hold on;

% Plot for PLO
semilogy(Convergence_curve_PLO_all, ':', 'LineWidth', 2.5, 'Color', [0.4660, 0.6740, 0.1880]); % Green dotted with triangle-down markers
hold on;

% Plot for PO
semilogy(Convergence_curve_PO_all, '-', 'LineWidth', 2.5, 'Color', [0.3010, 0.7450, 0.9330]); % Light blue solid with pentagram markers
hold on;

% Plot for ETO
semilogy(Convergence_curve_ETO_all, '--', 'LineWidth', 2.5, 'Color', [0.6350, 0.0780, 0.1840]); % Dark red dashed with hexagram markers
hold on;

% Plot for ESMA
semilogy(Convergence_curve_ESMA_all, '-.', 'LineWidth', 2.5, 'Color', [0.4940, 0.1840, 0.5560]); % Purple solid with star markers
hold on;

% Plot for AOSMA
semilogy(Convergence_curve_AOSMA_all, '-', 'LineWidth', 2.5, 'Color', [0.4660, 0.6740, 0.1880]); % Green solid with 'x' markers
hold on;

% Plot for LSMA
semilogy(Convergence_curve_LSMA_all, '-.', 'LineWidth', 2.5, 'Color', [0.8500, 0.3250, 0.0980]); % Red dotted with '+' markers
hold on;


xlabel('Iterations');
ylabel('Best score');


legend( 'ASMA','AO','SMA','QIO','PLO','PO','ETO','ESMA','AOSMA', 'LSMA');
title(['F1']);
box on; axis tight; grid on;
current_ylim = ylim;
ylim([current_ylim(1)*0.85, current_ylim(2)]);
% Save the figure
saveas(gcf, 'convergence_curves_ExampleF1_F23.png');
