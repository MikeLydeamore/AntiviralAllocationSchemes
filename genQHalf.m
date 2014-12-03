function [Q, stateList]=genQHalf(var,k)
%Generates the Q matrix for a single household of size k without antiviral
%intervention

%Generate state list
counter=1;
for s=1:k+1
	for e=1:k+1
		for i=1:k+1
			if (s-1)+(e-1)+(i-1)<=k
				stateList(:,counter)=[s-1;e-1;i-1];
				counter=counter+1;
			end
		end
	end
end
[m, noStates]=size(stateList);
Q=sparse(noStates,noStates);

%Infection events (S>0, I>0)
infectionEvents=intersect(find(stateList(1,:)>0),find(stateList(3,:)>0));

%Move to state (S-1, E+1, I) at rate (1-rho)(1-tau)beta*S*I/(k-1)
for i=1:length(infectionEvents)
	currState=infectionEvents(i);
	sminusone = find(stateList(1,:)==stateList(1,currState)-1);
	eplusone = find(stateList(2,:)==stateList(2,currState)+1);
	currenti = find(stateList(3,:)==stateList(3,currState));
	moveTo=intersect(intersect(sminusone,eplusone),currenti);
    if k~=1
        Q(currState,moveTo)=(1-var.rho)*(1-var.tau)*var.beta*stateList(1,currState)*stateList(3,currState)/(k-1);
    else
        Q(currState,moveTo)=(1-var.rho)*(1-var.tau)*var.beta*stateList(1,currState)*stateList(3,currState)/(k);
    end;
end

%Shedding events (E>0)
sheddingEvents=find(stateList(2,:)>0);

%Move to state (S, E-1, I+1) at rate sigma*E
for i=1:length(sheddingEvents)
	currState=sheddingEvents(i);
	currents = find(stateList(1,:)==stateList(1,currState));
	eminusone = find(stateList(2,:)==stateList(2,currState)-1);
	iplusone = find(stateList(3,:)==stateList(3,currState)+1);
	moveTo=intersect(intersect(currents,eminusone),iplusone);
	Q(currState,moveTo)=var.sigma*stateList(2,currState);
end

%Recovery events (I>0)
recoveryEvents=find(stateList(3,:)>0);

%Move to state (S, E, I-1) at rate (1+eta)gamma*I
for i=1:length(recoveryEvents)
	currState=recoveryEvents(i);
	currents = find(stateList(1,:)==stateList(1,currState));
	currente = find(stateList(2,:)==stateList(2,currState));
	iminusone = find(stateList(3,:)==stateList(3,currState)-1);
	moveTo=intersect(intersect(currents,currente),iminusone);
	Q(currState,moveTo)=(1+var.eta).*var.gamma.*stateList(3,currState);
end

%Leave the diagonals until the last step.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Generate Q                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add in diagonals
rowSum=sum(Q,2);
[m, noStates]=size(Q);
for i=1:noStates
	Q(i,i)=-rowSum(i);
end

end