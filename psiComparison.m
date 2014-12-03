psiVect=0:0.0002:1;

FES=zeros(1,length(psiVect));
peakTime=zeros(1,length(psiVect));
%peakSize=zeros(1,length(psiVect));
growth=zeros(1,length(psiVect));

stateList=[];
houseSizes=[];
for k=1:length(var.pi_k)
    [Q, s]=genQ(var,k,0);
    stateList=[stateList, s];
    houseSizes=[houseSizes, k*ones(1,length(s))];
end

for i=1:length(psiVect)
	fprintf('Current psi: %d\n',psiVect(i));
	var.psi=0;
	ic=init_conds_mh(var);
	var.psi=psiVect(i);
	
	[tout, yout]=selfConstMH(var,ic,initialCondition);
	
	FES(i) = stateList(3,:)*yout(end,:)'*var.N;
	
	infAtT=stateList(3,:)*yout'*var.N;
	[peakSize(i), index]=max(infAtT);
	peakTime(i)=tout(index);
end