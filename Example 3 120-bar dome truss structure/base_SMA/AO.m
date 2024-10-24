
%--------------------------------------------------------------------------
% 🌿🦟 Artemisinin Optimization (AO) 🦟🌿
% Artemisinin Optimization based on Malaria Therapy: Algorithm and Applications to Medical Image Segmentation
% 
% 🎨 Authors: Chong Yuan 📧 yc18338414794@163.com, Dong Zhao,
%            Ali Asghar Heidari 📧 as_heidari@ut.ac.ir, aliasghar68@gmail.com,
%            Lei Liu, Yi Chen, Zongda Wu, Huiling Chen 📧 chenhuiling.jlu@gmail.com
%
% 📅 Last update: May 1, 2024

% 🌟 Artemisinin Optimization (AO) source codes (version 1.0) 🌟
%
% Developed in MATLAB R2023a
%
% 📖 Reference: "Artemisinin Optimization based on Malaria Therapy: Algorithm and Applications to Medical Image Segmentation", Displays, Elsevier, 2024
%--------------------------------------------------------------------------

function [bestfitness,Leader_pos,Convergence_curve]=AO(N,MaxFEs,lb,ub,dim,fobj) %ASMA(N, 12500, lb, ub, dim, fobj);
% Initialization parameters
FEs=0;
it=1;

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
%% Main loop
while FEs<=MaxFEs
    tic; % Start timer for iteration  
    K= 1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    E =1*exp(-4*(FEs/MaxFEs));
    %
    for i=1: N
        Fitnorm(i)= (Fitness(i)-min(Fitness))/(max(Fitness)-min(Fitness));
        for j=1:dim
            if rand<K
                if rand<0.5
                    New_pop(i,j) = pop(i,j)+E.*pop(i,j)*(-1)^FEs;
                else
                    New_pop(i,j) = pop(i,j)+E.*best(j)*(-1)^FEs;
                end
            else
                New_pop(i,j)=pop(i,j);
            end
            if rand<Fitnorm(i)
                A=randperm(N);
                beta=(rand/2)+0.1;
                New_pop(i,j)=pop(A(3),j)+beta.*(pop(A(1),j)-pop(A(2),j)); 
            end
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
   % disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(Convergence_curve(it)) ': Time =' num2str(time_per_iteration(it)) ]);
 
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

function z=Transborder_reset(z,ub,lb,dim,best)
for j=1:dim
    if z(j)>ub || z(j)<lb
        
        z(j)=best(j);
        
    end
end
end


% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

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
