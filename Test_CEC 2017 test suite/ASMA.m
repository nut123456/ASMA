%Artemisinin slime mold algorithm (ASMA)
% Author and programmer: 
% Mr. Arnut Sutha
% Center of Excellence in Applied Mechanics and Structures, Department of Civil Engineering, Chulalongkorn University, 10330 Bangkok, Thailand
% e-mail:       mynut2009@gmail.com
% Researchgate: https://www.researchgate.net/profile/Arnut_Sutha
%_____________________________________________________________________________________________________    
function [bestfitness,Leader_pos,Convergence_curve]=ASMA(N,MaxFEs,lb,ub,dim,fobj) %SMA(N,Max_iter,lb,ub,dim,fobj)
% Initialization parameters
FEs=0;
it=1;
tFitness = inf*ones(1,N);%record the fitness of all slime mold
%% Initialization of the solution set
pop=initialization(N,dim,ub,lb);
%Calculate the fitness value of the initial solution set
for i=1:N
    Fitness(i)=fobj(pop(i,:));
    FEs=FEs+1;
end
[fmin,x]=min(Fitness);

%Container
New_pop=zeros(N,dim);
Fitnorm=zeros(1,N);
Convergence_curve=[];
%Record the current optimal solution
best=pop(x,:);
bestfitness=fmin;
time_per_iteration = zeros(1, MaxFEs); % Preallocate array for time per iteration
delta_threshold = 1e-4; % δ^t
%% Main loop
while FEs<=MaxFEs
    tic; % Start timer for iteration  
    [SmellOrder,SmellIndex] = sort(Fitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
    %calculate the fitness weight of each slime mold
    for i=1:N
        for j=1:dim
            if i<=(N/2)  %Eq.(2.5)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestfitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestfitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    E=1*exp(-4*(FEs/MaxFEs)^(1/2));
    a = 1;
    b = 0.003;
    % c = 0.05;
    d = 0.01; %0.3;
     % IMF inertia weight calculation
    % t = FEs/MaxFEs; % normalized iteration number (0 to 1)
    %E1 = a * exp(-b * FEs) + d; % IMF inertia weight equation
    vb = unifrnd(-E,E,1,dim);  %Eq.(2.3)
    Fitness_threshold = mean(Fitness);  % ใช้ค่าเฉลี่ยของ Fitness เป็นเกณฑ์
    for i=1: N
        Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));
        for j=1:dim
            if rand<K
                if rand<0.5
                    New_pop(i,j) = pop(j,j) + E .* sin((2*pi)*rand()) * pop(i,j) * (-1)^FEs;
                else
                    New_pop(i,j) = pop(i,j) + E .* cos((2*pi)*rand()) * best(j) * (-1)^FEs;
                end
            else
               New_pop(i,j)=pop(i,j);
                % New_pop(i,j) = (ub-lb)*rand+lb;
            end
            if rand<Fitnorm(i)                
                AA=randperm(N);   %randperm(N) สุ่มเรียงลำดับตัวเลขทั้งหมดจาก 1 ถึง N โดยไม่ซ้ำ.
                New_pop(i,j)=pop(AA(3),j)+ vb(j)*(weight(i,j)*pop(AA(1),j)-pop(AA(2),j));
            else
                    New_pop(i,j) = vb(j)*pop(i,j);
            end
        end
      
        % Exponential Decrease Adaptive Mutation
        P_max = 0.5; % Initial mutation probability
        P_min = 0.1; % Minimum mutation probability
        tau = MaxFEs / log(P_max / P_min); % Time constant for exponential decrease
        % Exponential Decrease Adaptive Mutation
        mutation_prob = P_max * exp(-FEs / tau);% exp(-FEs/ tau)% Exponential decrease formula
        if rand() < mutation_prob
            mutation_index = randi([1 N]);
            New_pop(mutation_index, :) = lb + (ub - lb) .* rand(1, dim);
        end
        New_pop(i,:)=Mutation(New_pop(i,:),pop(i,:),best,dim);
        New_pop(i,:)=Transborder_reset(New_pop(i,:),ub,lb,dim,best);
        
        tFitness=fobj(New_pop(i,:));
        FEs=FEs+1;
        if tFitness<Fitness(i)
            pop(i,:)= New_pop(i,:);
            Fitness(i)=tFitness;
        end
    end
    [fmin,x]=min(Fitness);
    if fmin<bestfitness
        best=pop(x,:);
        bestfitness=fmin;
    end
    Convergence_curve(it)=bestfitness;
    Leader_pos=best;
    time_per_iteration(it) = toc; % End timer and record elapsed time
    disp(['Iteration ' num2str(it) ':   Best Cost(ASMA) = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ]);
    
    %% Relative difference convergence check
    P_best = min(Fitness);        % Best penalized fitness
    P_mean = mean(Fitness);       % Mean penalized fitness
    delta = abs((P_mean - P_best)/P_best);  % Eq. (21)
    if delta < delta_threshold
        disp(['Terminated due to convergence: delta = ' num2str(delta) ', threshold = ' num2str(delta_threshold)]);
        break;
    end

    it=it+1;
end
disp(['Time(s) = ' num2str(sum(time_per_iteration))]);
end

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

% 
% function z=Transborder_reset(z,ub,lb,dim,best)
% for j=1:dim
%     if z(j)>ub || z(j)<lb
% 
%         z(j)=best(j);
% 
%     end
% end
% end

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
% function Positions = initialization(SearchAgents_no,dim,ub,lb)
% % This function initialize the first population of search agents
% Boundary_no = size(ub,2); % Numnber of boundaries
% % If the boundaries of all variables are equal and user enter a signle
% % Number for both ub and lb
% if Boundary_no == 1
%     Positions = rand(SearchAgents_no,dim).*(ub-lb)+lb;
% end
% % If each variable has a different lb and ub
% if Boundary_no > 1
%     for i = 1:dim
%         ub_i = ub(i);
%         lb_i = lb(i);
%         Positions(:,i) = rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
%     end
% end
% end

