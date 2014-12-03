function [infected, eventTime, noAV, totalInfected, popSize, houseInfo]=sim_run_expexp(var,initCond,stateList,houseSizes)
%Runs a simulation for an epidemic for an exp-exp type problem.
%Kappa and zeta are the time, not the rates.

%% Setup
%Distribute household sizes (must add up to N)
housePop=zeros(1,length(var.h_k));
for i=1:length(var.h_k)-1
    housePop(i)=round(var.h_k(i)*var.N);
end
housePop(end)=var.N-sum(housePop);
popSize=housePop*(1:length(housePop))';

%Set up the initial discrete conditions
housesEachType=discreteConditions(var,initCond,stateList);
nHouses=sum(housesEachType);
housesEachType=cumsum([0 housesEachType]);

%Maximum iterations is now 4*popSize+N, as there is a transition in each
%household where they get the antivirals.
maxIter=4*popSize+var.N;

%Load in the population
popStatus=zeros(nHouses,5); %S-E-I-R + Antiviral Status

cumHouseSizes=[0 cumsum(housePop)];
for i=2:length(cumHouseSizes)
    popStatus(cumHouseSizes(i-1)+1:cumHouseSizes(i),1)=i-1;
end

%Set up the initCond
houseSizeDivision=zeros(housesEachType(end),1);
for i=1:length(housesEachType)-1
    popStatus(housesEachType(i)+1:housesEachType(i+1),1)=stateList(1,i);
    popStatus(housesEachType(i)+1:housesEachType(i+1),2)=stateList(2,i);
    popStatus(housesEachType(i)+1:housesEachType(i+1),3)=stateList(3,i);
    popStatus(housesEachType(i)+1:housesEachType(i+1),4)=houseSizes(i)-(stateList(1,i)+stateList(2,i)+stateList(3,i));
	 popStatus(housesEachType(i)+1:housesEachType(i+1),5)=stateList(4,i);
     if houseSizes(i)==1
         houseSizeDivision(housesEachType(i)+1:housesEachType(i+1))=1;
     else
         houseSizeDivision(housesEachType(i)+1:housesEachType(i+1))=houseSizes(i)-1;
     end
end

%Preallocation is taken care of in the initial condition now

%Preallocate the antivirals
%phi_k(k) is the proportion of households of size k that are preallocated AVs.
%preAllocation=round(var.phi_k.*housePop);
%Load in antivirals
%cumSizes=[0 cumsum(housePop)];
%for i=1:length(cumSizes)-1
%    popStatus(cumSizes(i)+1:cumSizes(i)+preAllocation(i),5)=1;
%end

eventTime=zeros(1,maxIter);

infected=zeros(1,maxIter);
noAV=sum((popStatus(:,5)~=0).*sum(popStatus(:,1:4),2).*var.dosage);
AVTime=zeros(var.N,1);
houseInfo=zeros(maxIter,5);

maxT=1000;
counter=1;
t=0;


%% Simulation
%Rate vector goes [Infection ... Progression ... Recovery ... A.V. Intro A.V. Removal]
rates=zeros(6*var.N,1);
while t<maxT && sum(popStatus(:,2)+popStatus(:,3))>0
    infected(counter)=sum(popStatus(:,3));
    if infected(counter)>popSize
        error('Oh dear...');
	 end

	 if mod(counter,500)==0
        clc;
        fprintf('Total infected: %d\n',popSize-sum(popStatus(:,1)));
        fprintf('Current Infected: %d\n',infected(counter));
		  fprintf('AV''s allocated: %d\n',noAV);
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Rates                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %hasAVVector=popStatus(:,5)==1;
	 reduce=popStatus(:,5)==1;
	 hasAVVector=reduce;
	 %fullReduction=find(((popStatus(:,3)+popStatus(:,4))~=0)&popStatus(:,5)==1); %Full on reductions
	 %earlyReduction=find(((popStatus(:,3)+popStatus(:,4))==0)&popStatus(:,5)==1); %Early starters
    %Infection - Internal
    %If household has antivirals, then beta will be modified by
    %(1-tau)(1-rho)
    modifiers = (1-hasAVVector.*var.tau).*(1-hasAVVector.*var.rho);
    
    rates(1:var.N)=var.beta./(houseSizeDivision).*popStatus(:,1).*popStatus(:,3).*modifiers;
    
    %Infection - External
    %If the external household has antivirals, then we have (1-tau)
    %If the internal household has antivirals, then we have (1-rho)
    %Need the # of infected from houses with antirivals and without
    %Then the external force is alpha*(1-tau)*I_a + alpha*I_{no a}
    noInfWithAV=(hasAVVector)'*popStatus(:,3);
    noInfWithoutAV=not(hasAVVector)'*popStatus(:,3);
    
    tfi=var.alpha*((1-var.tau)*noInfWithAV+noInfWithoutAV)./popSize;
    rates(1:var.N)=rates(1:var.N)+tfi*popStatus(:,1).*(1-hasAVVector.*var.rho);
    
    %Latent Progression
    %Progression always at the same rate, sigma*E
    rates(var.N+1:2*var.N)=var.sigma*popStatus(:,2);
    
    %Recovery
    %Recovery is at rate gamma*I if there's no antivirals, and
    %(1+ita)*gamma*I if there is antivirals
    rates(2*var.N+1:3*var.N)=(1+hasAVVector*var.eta).*var.gamma.*popStatus(:,3);
    
    %Antiviral introduction
    %Constant rate 1/zeta into households who have at least 1 infection event
    %and don't already have antivirals
    rates(3*var.N+1:4*var.N)=1/var.zeta*((popStatus(:,3)+popStatus(:,4))>0).*(popStatus(:,5)==0);
    
    %Antiviral removal
    %Constant rate 1/kappa for households who already have the antivirals
    rates(4*var.N+1:5*var.N)=1/var.kappa*(hasAVVector);
	 
	 %False starts
	 rates(5*var.N+1:6*var.N)=(popStatus(:,5)==4).*var.psi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         Choosing Events                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumRates=sum(rates);
    %Normalise
    normRateVector=rates./sumRates;
    
    %Choose event
    u=rand();
    eventOccured=find(u<cumsum(normRateVector),1);
    
    %Calculate event time
    eventTime(counter+1)=exprnd(1/sumRates);
    t=t+eventTime(counter+1);
    
    %Determine event and household
    household=mod(eventOccured,var.N);
    if household==0
        household=var.N;
    end
    if eventOccured<var.N+1 %Infection
        popStatus(household,1)=popStatus(household,1)-1; %S -> S-1
        popStatus(household,2)=popStatus(household,2)+1; %E -> E+1
        
    elseif eventOccured < 2*var.N+1 %Progression
        popStatus(household,2)=popStatus(household,2)-1; %E -> E-1
        popStatus(household,3)=popStatus(household,3)+1; %I -> I+1
			if popStatus(household,5)==4 %Preallocated
				popStatus(household,5)=1;
				noAV=noAV+var.dosage*sum(popStatus(household,1:4));
				AVTime(household)=t;
			end
        
    elseif eventOccured < 3*var.N+1 %Recovery
        popStatus(household,3)=popStatus(household,3)-1; %I -> I-1
        popStatus(household,4)=popStatus(household,4)+1; %R -> R+1
        
    elseif eventOccured < 4*var.N+1 %Antivirals introduced

        popStatus(household,5)=1;
		  noAV=noAV+var.dosage*sum(popStatus(household,1:4));
        AVTime(household)=t;
		  
		  if (var.maxAV-noAV)<var.k
			  var.zeta=Inf;
		  end
        
	 elseif eventOccured < 5*var.N+1 %Antivirals finished
		 popStatus(household,5)=2;
		 
	 else %False start
		 popStatus(household,5)=1;
		 noAV=noAV+var.dosage*sum(popStatus(household,1:4));
		 AVTime(household)=t;
		 %fprintf('False start in household %d\n',household);
    end
    houseInfo(counter,:)=popStatus(500,:);
    counter=counter+1;
end
infected(counter)=sum(popStatus(:,3));
eventTime=cumsum(eventTime);
totalInfected=popSize-sum(sum(popStatus(:,4)));
end
