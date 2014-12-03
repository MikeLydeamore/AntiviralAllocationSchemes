function dPdt=selfConstDPMH(var,Q,Q2,QWithoutAV,QVaccinationOnly,stateList,reduceStates,houseSizes,p,t)
%Function which performs time dependent calculations for Self-Consistent
%Field Method.
%
%Syntax: dPdt=selfConstDPMH(var,Q,Q2,QWithoutAV,stateList,reduceStates,houseSize,p,t)
%
%Designed for use strictly with selfConstMH
%
%Michael Lydeamore - The University of Adelaide - 2014

%Feed in almost everything, then just modify by the current proportion
%infected.

%Find the mean household size:
meansize=1:length(var.pi_k);
meansize=meansize*var.h_k';
%Find the external infection rate
proportion_infected=p'* ((ones(1,length(stateList)) .* (1-((stateList(4,:)==1).*var.tau))) .*stateList(3,:))'/meansize;
externalInfectionMatrix=var.alpha.*proportion_infected.*Q2;


%Reduce for antivirals
externalInfectionMatrix(reduceStates,:)=(1-var.rho)*externalInfectionMatrix(reduceStates,:);
%%%

%Get current number of Antivirals
haveAV = (stateList(4,:)==1)|(stateList(4,:)==2)|(stateList(4,:)==4)|(stateList(4,:)==6);
currentAV=(haveAV*(houseSizes.*var.N.*p')');


if (var.maxAV-currentAV) < var.k
    fullQ=QWithoutAV+externalInfectionMatrix;
else
    fullQ=Q+externalInfectionMatrix;
end

if (t>var.vaccTime)
    fullQ=fullQ+QVaccinationOnly;
end

dPdt=p'*fullQ;
dPdt=dPdt';

end