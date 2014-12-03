function [tout yout]=selfConstMH(var,initialCondition,tspan)
%Set's up and solves the Self-Consistent Field Method equations for a
%pandemic outbreak
%
%Syntax: [tout, yout]=selfConstMH(var,initialCondition,tspan)
%
%Dependencies: genQ
%              selfConstDPMH
%Does not strictly depend on init_conds_mh, however the initialCondition
%output from init_conds_mh is of the correct form for this function.
%
%Michael Lydeamore - The University of Adelaide - 2014

%variables;
var.zeta=1/var.zeta;
var.kappa=1/var.kappa;
Q=[];
megaStates=[];
statesLength=[];
QWithoutAV=[];
houseSizes=[];
QVaccinationOnly=[];
tempnu = var.nu;
for k=1:length(var.pi_k)
    
    var.nu=0; %Ignore vaccination
    
    %Construct regular household Q
    [Qh, stateList]=genQ(var,k,0);
    [mh, nh]=size(Qh);
    [mq, nq]=size(Q);
    topPadding=zeros(mq,nh);
    sidePadding=zeros(mh,nq);
    Q=[Q, topPadding ; sidePadding, Qh];
    statesLength=[statesLength, length(stateList)];
    megaStates=[megaStates, stateList];
    houseSizes=[houseSizes, k*ones(1,length(stateList))];
    
    %Generate household P without the ability to have antivirals.
    noAVVar=var;
    noAVVar.zeta=0;
    QWithoutAVH=genQ(noAVVar,k,0);
    [mh, nh]=size(QWithoutAVH);
    [mq, nq]=size(QWithoutAV);
    topPadding=zeros(mq,nh);
    sidePadding=zeros(mh,nq);
    QWithoutAV=[QWithoutAV, topPadding ; sidePadding, QWithoutAVH];
    
    %Vaccination only
    var.nu=tempnu;
    QVacc = genQ(var,k,0);
    setToZero = (stateList(4,:)~=5)&(stateList(4,:)~=6);
    QVacc(setToZero,setToZero)=0;
    [mh, nh]=size(QVacc);
    [mq, nq]=size(QVaccinationOnly);
    topPadding=zeros(mq,nh);
    sidePadding=zeros(mh,nq);
    QVaccinationOnly=[QVaccinationOnly, topPadding; sidePadding, QVacc];
    
end
statesLength=[1, cumsum(statesLength)];
%Construct the external force matrix
[m, n]=size(Q);
Q2=sparse(m,n);
%There is a possible external infection for every state except those with 0
%susceptibles.
externalInfectionStates=find(megaStates(1,:)~=0);
for l=1:length(externalInfectionStates)
    hsize=find(externalInfectionStates(l)<=statesLength,1);
    s=megaStates(1,externalInfectionStates(l));
    e=megaStates(2,externalInfectionStates(l));
    i=megaStates(3,externalInfectionStates(l));
    a=megaStates(4,externalInfectionStates(l));

	 transitionTo=find(megaStates(1,:)==s-1&megaStates(2,:)==e+1&megaStates(3,:)==i&megaStates(4,:)==a);

    if length(transitionTo)>1
        index=find(transitionTo>statesLength(hsize-1)&transitionTo<=statesLength(hsize));
        transitionTo=transitionTo(index);
    end
    Q2(externalInfectionStates(l),transitionTo)=s;
end
%Balance out the row sums
for i=1:length(Q2)
    Q2(i,i)=-sum(Q2(i,:));
end

%Output from init_conds
P=initialCondition;

%Get states for antiviral reduction
reduceStates=[];
for k=1:length(var.pi_k)
    st=find(megaStates(4,:)==1);
    st=st(st>statesLength(k)&st<statesLength(k+1));
    reduceStates=[reduceStates, st];
end

    options=odeset('OutputFcn',@(t,y,flag) destatus(t,y,flag,var,megaStates));


%fix up the diagonals for the vaccination
rowSum=sum(QVaccinationOnly,2);
[m, noStates]=size(QVaccinationOnly);
for i=1:noStates
    QVaccinationOnly(i,i)=-rowSum(i);
end

dPdt=@(t,y) selfConstDPMH(var,Q,Q2,QWithoutAV,QVaccinationOnly,megaStates,reduceStates,houseSizes,y,t);

[tout, yout]=ode45(dPdt, tspan, P);

It=var.N*megaStates(3,:)*yout';

logInf=It;
logInf(logInf>0)=log(logInf(logInf>0));


%Total number of infected:
totalRecovered=yout(end,:)*var.N*(k-sum(megaStates(1:3,:))');

