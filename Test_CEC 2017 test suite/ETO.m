%%Exponential-Trigonometric Optimization Algorithm

function [Destination_fitness,Destination_position,Convergence_curve]=ETO(N,Max_Iter,LB,UB,Dim,Fobj)
Destination_position=zeros(1,Dim);
Destination_fitness=inf;
Destination_position_second=zeros(1,Dim);
Convergence_curve=zeros(1,Max_Iter);
Position_sort = zeros(N,Dim);
%Initialize ETO parameters
b=1.55;
CE=floor(1+(Max_Iter/b));
T=floor(1.2+Max_Iter/2.25);
CEi=0;
CEi_temp=0;
UB_2=UB;
LB_2=LB;
%Initialize the set of random solutions
X=initialization(N,Dim,UB,LB);
Objective_values = zeros(1,size(X,1));
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=Fobj(X(i,:));
    if Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
end
Convergence_curve(1)=Destination_fitness;
t=2; 
%Main loop
while t<=Max_Iter    
    for i=1:size(X,1) % in i-th solution
        
        for j=1:size(X,2) % in j-th dimension
           
            %update A by using Eq. (17)
            d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (1 - t / Max_Iter));
            d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (1 - t / Max_Iter));
          
            CM=(sqrt(t/Max_Iter)^tan(d1/(d2)))*rand()*0.01; 
            % enter the bounded search strategy
            if t==CEi
                UB_2=Destination_position(j)+(1-t/Max_Iter)*abs(rand()*Destination_position(j)-Destination_position_second(j))*rand();
                LB_2=Destination_position(j)-(1-t/Max_Iter)*abs(rand()*Destination_position(j)-Destination_position_second(j))*rand();
                if UB_2>UB
                    UB_2=UB;
                end
                if LB_2<LB
                    LB_2=LB;
                end
                X=initialization(N,Dim,UB_2,LB_2);                
                CEi_temp=CEi;
                CEi=0;
            end
        % the first phase of exploration and exploitation    
        if t<=T%3.6-3.62  
            q1=rand();
            q3=rand();
            q4=rand();
            if CM>1
                d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q1));
                d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q1));
                alpha_1=rand()*3*(t/Max_Iter-0.85)*exp(abs(d1/d2)-1);
                if q1<=0.5
                    X(i,j)=Destination_position(j)+rand()*alpha_1*abs(Destination_position(j)-X(i,j));
                else
                    X(i,j)=Destination_position(j)-rand()*alpha_1*abs(Destination_position(j)-X(i,j));  
                end                
            else
                d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q3));
                d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q3));
                alpha_3=rand()*3*(t/Max_Iter-0.85)*exp(abs(d1/d2)-1.3);
                if q3<=0.5
                    X(i,j)=Destination_position(j)+q4*alpha_3*abs(rand()*Destination_position(j)-X(i,j));
                else
                    X(i,j)=Destination_position(j)-q4*alpha_3*abs(rand()*Destination_position(j)-X(i,j));  
                end
            end
        else
            % the second phase of exploration and exploitation
            q2=rand();
            alpha_2=rand()*exp(tanh(1.5*(-t/Max_Iter-0.75) - rand()));
            if CM<1
                d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q2));
                d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_Iter * (q2));
                X(i,j)= X(i,j)+exp(tan(abs(d1/d2))*abs(rand()*alpha_2*Destination_position(j)-X(i,j)));
            else
                if q2<=0.5
                    X(i,j)=X(i,j)+3*(abs(rand()*alpha_2*Destination_position(j)-X(i,j)));
                else
                    X(i,j)=X(i,j)-3*(abs(rand()*alpha_2*Destination_position(j)-X(i,j)));  
                end
            end 
        end
        end
        CEi=CEi_temp;
    end
    for i=1:size(X,1)         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>UB_2;
        Flag4lb=X(i,:)<LB_2;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+(UB_2+LB_2)/2.*Flag4ub+LB_2.*Flag4lb;        
        % Calculate the objective values
        Objective_values(1,i)=Fobj(X(i,:));
%         % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    % find the second solution
    if t==CE
        CEi=CE+1;
        CE=CE+floor(2-t*2/(Max_Iter-CE*4.6)/1);
        temp = zeros(1,Dim);
        temp2 = zeros(N,Dim);
        %sorting
        for i=1:(size(X,1)-1)
	        for j=1:(size(X,1)-1-i)
		        if Objective_values(1,j) > Objective_values(1,j+1)
			        temp(1,j) = Objective_values(1,j);
			        Objective_values(1,j) = Objective_values(1,j+1);
			        Objective_values(1,j+1) = temp(1,j);
			        temp2(j,:) = Position_sort(j,:);
			        Position_sort(j,:) = Position_sort(j+1,:);
			        Position_sort(j+1,:) = temp2(j,:);   
                end
            end
        end
        Destination_position_second=Position_sort(2,:);%the second solution
    end
    Convergence_curve(t)=Destination_fitness;
    t=t+1;
end

function X=initialization(SearchAgents_no,Dim,UB,LB)
Boundary_no= size(UB,2); 
% If the boundaries of all variables are equal
if Boundary_no==1
    X=rand(SearchAgents_no,Dim).*(UB-LB)+LB;  
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:Dim
        UB_i=UB(i);
        LB_i=LB(i);
        X(:,i)=rand(SearchAgents_no,1).*(UB_i-LB_i)+LB_i;
    end
end