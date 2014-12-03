function [initialCondition, growth, Jacob]=init_conds_mh(var)
%Determines the 'correct' initial conditions for immediate takeoff, using
%eigenvalues of the Jacobian.
%
%Syntax: [initialCondition, growthRate, Jacobian]=init_conds_mh(var);
%
%Dependencies:
%   genQ - Does not require Q Matrix, but does require the list of states.
%   variables - Struct containing parameters for pandemic model.
%
%Michael Lydeamore - The University of Adelaide - 2014
%Edited 14/04/2014
%


%%% Disable this if Q matrix (and other things) is in the current folder
%addpath('..');

prop_inf=2*10^-3;

JacIndex=0;
Jacob=[];
indexes=[];
megaStates=[];
AVIndexes=[];
statesLength=[];
for k=1:length(var.pi_k)
    
    [Q, stateList]=genQ(var,k,0);
    
    full_s = find(stateList(1,:)==k&stateList(2,:)==0&stateList(3,:)==0&stateList(4,:)==0);
    full_s_withAV = find(stateList(1,:)==k&stateList(2,:)==0&stateList(3,:)==0&stateList(4,:)==4);
    
    
    %Build Jacobian:
    %Matrix of recovery events:
    Rec=sparse(zeros(length(stateList)));
    recoveryEvents=find(stateList(3,:)>0);
    
    for i=1:length(recoveryEvents)
        currState=recoveryEvents(i);
        moveTo= find(stateList(1,:)==stateList(1,currState)&...
            stateList(2,:)==stateList(2,currState)&...
            stateList(3,:)==stateList(3,currState)-1&...
            stateList(4,:)==stateList(4,currState));
        Rec(currState,moveTo)=1;
        
    end
    rowSum=sum(Rec,2);
    [m, noStates]=size(Rec);
    for i=1:noStates
        Rec(i,i)=-rowSum(i);
    end
    DRec=sparse(diag(var.gamma*stateList(3,:)));
    
    JacRec = Rec'*DRec';
    
    %Matrix of all infection events
    Infec=sparse(zeros(length(stateList)));
    infectionEvents=find(stateList(1,:)>0);
    
    for i=1:length(infectionEvents)
        currState=infectionEvents(i);
			  moveTo= find(stateList(1,:)==stateList(1,currState)-1&...
				  stateList(2,:)==stateList(2,currState)+1&...
				  stateList(3,:)==stateList(3,currState)&...
				  stateList(4,:)==stateList(4,currState));
        Infec(currState,moveTo)=1;
    end
    
    rowSum=sum(Infec,2);
    [m, noStates]=size(Infec);
    for i=1:noStates
        Infec(i,i)=-rowSum(i);
    end
    
    SVec = stateList(1,:).*stateList(3,:);
    modifier=(1-((stateList(4,:)==1).*var.tau)).*(1-(((stateList(4,:)==1)).*var.rho));
    
    if k==1
        DInInf = sparse(diag(modifier.*var.beta.*SVec));
    else
        DInInf = sparse(diag(modifier.*(var.beta/(k-1)).*SVec));
    end
    
    JacInInf = Infec'*DInInf';
    
    
    %Progression matrix: Infection events
    Prog = sparse(zeros(length(stateList)));
    progressionEvents=find(stateList(2,:)>0);
    
    for i=1:length(progressionEvents)
        currState=progressionEvents(i);
		  if stateList(4,currState)~=4
			  moveTo= stateList(1,:)==stateList(1,currState)&...
				  stateList(2,:)==stateList(2,currState)-1&...
				  stateList(3,:)==stateList(3,currState)+1&...
				  stateList(4,:)==stateList(4,currState);
		  else
			  moveTo= stateList(1,:)==stateList(1,currState)&...
				  stateList(2,:)==stateList(2,currState)-1&...
				  stateList(3,:)==stateList(3,currState)+1&...
				  stateList(4,:)==1;
		  end
        Prog(currState,moveTo)=1;
    end
    rowSum=sum(Prog,2);
    [m, noStates]=size(Prog);
    for i=1:noStates
        Prog(i,i)=-rowSum(i);
    end
    
    DProg = sparse(diag(var.sigma*stateList(2,:)));
    JacProg = Prog' * DProg';
    
    
    %%% Antivirals
    % For transitions into states with antivirals, we have an 'almost-diagonal'
    % structure. There is a -sigma in the 'top' section whenever there has been
    % at least one infection, and a sigma in the corresponding 'middle' state.
    transitionPossible=find((stateList(1,:)+stateList(2,:))~=k&stateList(4,:)==0);
    intoAV=sparse(length(stateList),length(stateList));
    d=sum(stateList(4,:)==4);
    shift=(length(stateList)-d-2)/3; %Minus 2 Vaccination states
    for i=1:length(transitionPossible)
        intoAV(transitionPossible(i),transitionPossible(i))=-1/var.zeta; %Flow out
        intoAV(transitionPossible(i)+shift, transitionPossible(i))=1/var.zeta; %Flow in
    end
    
    transitionPossible=find(stateList(4,:)==1);
    outOfAV=sparse(length(stateList),length(stateList));
    for i=1:length(transitionPossible)
        outOfAV(transitionPossible(i),transitionPossible(i))=-1/var.kappa; %Flow out
        outOfAV(transitionPossible(i)+shift,transitionPossible(i))=1/var.kappa; %Flow in
	 end
	 
	 
	 %False taking - move from (s,e,i,4) to (s,e,i,1)
     falseStates=find(stateList(4,:)==4);
	 falseTaking=sparse(length(stateList),length(stateList));
	 for i=1:length(falseStates)
         currState=falseStates(i);
		 moveTo=find(stateList(1,:)==stateList(1,currState)&...
			 stateList(2,:)==stateList(2,currState)&...
			 stateList(3,:)==stateList(3,currState)&...
			 stateList(4,:)==1);
		 falseTaking(moveTo,currState)=var.psi;
		 falseTaking(currState,currState)=falseTaking(currState,currState)-(var.psi);
	 end    
    
    
    JacobTemp = JacRec+JacInInf+JacProg    + intoAV + outOfAV + falseTaking; 

    
    %Pad out the matrix to be of the right dimension.
    [mj, nj]=size(Jacob);
    [mt, nt]=size(JacobTemp);
    paddingAbove=sparse(zeros(mj,nt));
    paddingBeside=sparse(zeros(mt,nj));
    Jacob=[Jacob, paddingAbove; paddingBeside, JacobTemp];
    indexes = [indexes, length(megaStates)+find(stateList(1,:)==k&stateList(2,:)==0&stateList(3,:)==0&stateList(4,:)==0)];
    AVIndexes = [AVIndexes, length(megaStates)+find(stateList(1,:)==k&stateList(2,:)==0&stateList(3,:)==0&stateList(4,:)==4)];
    megaStates=[megaStates, stateList];
    statesLength=[statesLength, length(stateList)];
end
statesLength=[1, cumsum(statesLength)];

%External infection events:
DExInf = sparse(zeros(length(megaStates)));


for i=1:length(indexes)
    DExInf(indexes(i),:)=var.pi_k(i).*(1-var.phi_k(i)).*ones(1,length(megaStates)).*(1-(megaStates(4,:)==1).*var.tau).*var.alpha.*megaStates(3,:);
    DExInf(AVIndexes(i),:)=var.pi_k(i).*(var.phi_k(i)).*ones(1,length(megaStates)).*(1-(megaStates(4,:)==1).*var.tau).*var.alpha.*megaStates(3,:);
end

%Massive matrix of infection events
Infec=zeros(length(megaStates));
infectionEvents=find(megaStates(1,:)>0);
for ii=1:length(infectionEvents)
    hsize=find(infectionEvents(ii)<=statesLength,1);
    s=megaStates(1,infectionEvents(ii));
    e=megaStates(2,infectionEvents(ii));
    i=megaStates(3,infectionEvents(ii));
    a=megaStates(4,infectionEvents(ii));
    moveTo=find(megaStates(1,:)==s-1&megaStates(2,:)==e+1&megaStates(3,:)==i&megaStates(4,:)==a);
    if length(moveTo)>1
        index= moveTo>statesLength(hsize-1)&moveTo<=statesLength(hsize);
        moveTo=moveTo(index);
    end
    Infec(infectionEvents(ii),moveTo)=1;
end
rowSum=sum(Infec,2);
[m, noStates]=size(Infec);
for i=1:noStates
    Infec(i,i)=-rowSum(i);
end
JacExInf = Infec'*DExInf;

Jacob = Jacob+JacExInf;

opts.disp=0;
[v, d]=eigs(Jacob,1,'LR',opts);

eval=find(diag(d)>0);
growth=d(eval,eval);

%We have I(0) = S + eps*EVect, and I(0)*megaStates(3,:)=prop_inf, so
%rearranging...
eps = prop_inf/(v(:,eval)'*megaStates(3,:)');

initialCondition=zeros(1,length(megaStates));

initialCondition(indexes)=var.h_k.*(1-var.phi_k);
initialCondition(AVIndexes)=var.h_k.*var.phi_k;
initialCondition=initialCondition+eps*v(:,eval)';