%Artemisinin slime mold algorithm (ASMA)
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________           

function [bestfitness,Leader_pos,Convergence_curve]=ASMA(N,MaxFEs,lb,ub,dim,fobj) 
% Initialization parameters
FEs=0;
it=1;
tFitness = inf*ones(1,N);%record the fitness of all Artesinin slime mold
%% Initialization of the solution set
AS=initialization(N,dim,ub,lb);
%Calculate the fitness value of the initial solution set
for i=1:N
    Fitness(i)=fobj(AS(i,:));
    FEs=FEs+1;
end
[fmin,x]=min(Fitness);
%Container
New_AS=zeros(N,dim);
Fitnorm=zeros(1,N);
Convergence_curve=[];
%Record the current optimal solution
best=AS(x,:);
bestfitness=fmin;
time_per_iteration = zeros(1, MaxFEs); % Preallocate array for time per iteration
%% Main loop
while FEs<=MaxFEs
    tic; % Start timer for iteration  
    [SmellOrder,SmellIndex] = sort(Fitness);  %Eq.(8)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(8)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestfitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestfitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));  %Eq.(9)
    E=1*exp(-4*(FEs/MaxFEs));  %Eq.(12)
    vb = unifrnd(-E,E,1,dim);  %Eq.(12)
    for i=1: N
        Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));% Eq. (14)
        for j=1:dim
            if rand<K
                if rand<0.5
                    New_AS(i,j) = AS(i,j)+E.*AS(i,j)*(-1)^FEs; % Eq. (16)
                else
                    New_AS(i,j) = AS(i,j)+E.*best(j)*(-1)^FEs; % Eq. (16)
                end
            else
                New_AS(i,j)=AS(i,j);
            end
            if rand<Fitnorm(i)                
                AA=randperm(N);  %randperm(N) generates a random permutation of the integers from 1 to N without repetition
                New_AS(i,j)=AS(AA(3),j)+ vb(j)*(weight(i,j)*AS(AA(1),j)-AS(AA(2),j));% Eq. (17)
            end
        end
        % Exponential Decrease Adaptive Mutation %Eq.(18)
        P_max = 0.5; % Initial mutation probability
        P_min = 0.1; % Minimum mutation probability
        tau = MaxFEs / log(P_max / P_min); % Time constant for exponential decrease
        % Exponential Decrease Adaptive Mutation
        mutation_prob = P_max * exp(-FEs / tau); % Exponential decrease formula
        if rand() < mutation_prob
            mutation_index = randi([1 N]);
            New_AS(mutation_index, :) = lb + (ub - lb) .* rand(1, dim);
        end
        New_AS(i,:)=Mutation(New_AS(i,:),AS(i,:),best,dim);  %Eq.(15)
        New_AS(i,:)=Transborder_reset(New_AS(i,:),ub,lb,dim,best);
        
        tFitness=fobj(New_AS(i,:));
        FEs=FEs+1;
        if tFitness<Fitness(i)
            AS(i,:)= New_AS(i,:);
            Fitness(i)=tFitness;
        end
    end
    [fmin,x]=min(Fitness);
    if fmin<bestfitness
        best=AS(x,:);
        bestfitness=fmin;
    end
    Convergence_curve(it)=bestfitness;
    Leader_pos=best;
    time_per_iteration(it) = toc; % End timer and record elapsed time
    disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ]);
    it=it+1;
end
disp(['Time(s) = ' num2str(sum(time_per_iteration))]);
end

%% Mutation strategy
function z=Mutation(z,x,b,dim)
for j=1:dim
    if rand<0.05
        z(j)=x(j);
    end
    if rand<0.2
        z(j)=b(j);
    end
end
end

function z=Transborder_reset(z,ub,lb,dim,best)
    for j=1:dim
        if z(j) > ub || z(j) < lb
            z(j) = best(j);
        end
    end
end

%% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)
%Boundary_no = numel(ub); % number of boundaries
Boundary_no= size(ub,2); % numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end
