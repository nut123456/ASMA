%% The Adative Opposition Slime Mould Algorithm (AOSMA) function
function [Destination_fitness,bestPositions,Convergence_curve]=AOSMA(N,Max_iter,lb,ub,dim,fobj)

bestPositions=zeros(1,dim);
Destination_fitness=inf;%change this to -inf for maximization problems
AllFitness = inf*ones(N,1);%record the fitness of all slime mold
weight = ones(N,dim);%fitness weight of each slime mold
%Initialize the set of random solutions

X=initialization(N,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);

it=1;  %Number of iterations
lb=ones(1,dim).*lb; % lower boundary 
ub=ones(1,dim).*ub; % upper boundary
del=0.03; % parameter
time_per_iteration = zeros(1, Max_iter); % Preallocate array for time per iteration
for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));       
end


% Main loop
while  it <= Max_iter
    tic; % Start timer for iteration  
    oldfitness=AllFitness;
    [SmellOrder,SmellIndex] = sort(oldfitness);  %Eq.(7)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);

    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
  
    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(6)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    a = atanh(-(it/Max_iter)+1);   %Eq.(11)
    b = 1-it/Max_iter;             %Eq.(12)
    % Update the Position of search agents
    for i=1:N       %Eq.(13)
        if rand<del     
            X(i,:) = (ub-lb)*rand+lb;   %Eq.(13c)
        else
            p =tanh(abs(oldfitness(i)-Destination_fitness));  %Eq.(4)
            vb = unifrnd(-a,a,1,dim);  
            vc = unifrnd(-b,b,1,dim);
            A = randi([1,N]);  % one positions randomly selected from population
            r2 = rand();
            for j=1:dim
                if r2<p    %Eq.(13a)
                    X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*bestPositions(j)-X(A,j));
                else       %Eq.(13b)
                    X(i,j) = vc(j)*X(i,j);
                end
            end
        end
    end
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
    end
    for i=1:N
        if AllFitness(i)>oldfitness(i)              %Eq.(16)
            D(i,:)=min(X(i,:))+max(X(i,:))-X(i,:);  %Eq.(14)
            Flag4ub=D(i,:)>ub;
            Flag4lb=D(i,:)<lb;
            D(i,:)=(D(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            if fobj(D(i,:))<AllFitness(i)           %Eq.(15)
                X(i,:)=D(i,:);
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
    
    time_per_iteration(it) = toc; % End timer and record elapsed time
    %disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(Destination_fitness) ': Time =' num2str(time_per_iteration(it)) ]);
    it=it+1;
end
disp(['Time(s) = ' num2str(sum(time_per_iteration))]);

% Plot time per iteration
% figure(1);
% plot(1:2:Max_iter, time_per_iteration(1:2:Max_iter), '-or', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Time (s)');
% title('\fontsize{12}\bf Computation Time per Iteration');
% legend('\fontsize{10}\bf Time per Iteration');
% % Plot convergence curve (Best Cost over Iterations)
% figure(1);
% plot(1:2:Max_iter, Convergence_curve(1:2:Max_iter), '-ob', 'LineWidth', 1.4, 'MarkerSize', 3);
% xlabel('\fontsize{12}\bf Iteration');
% ylabel('\fontsize{12}\bf Best Cost');
% title('\fontsize{12}\bf Convergence Curve');
% legend('\fontsize{10}\bf Best Cost');


end