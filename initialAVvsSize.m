
%variables;
incVect=0:2000:100000;
delayVect=0.01:0.1:4;
preallocavs=cell(length(incVect),length(delayVect));
pret=cell(length(incVect),length(delayVect));
dynallocavs=cell(length(incVect),length(delayVect));
dynt=cell(length(incVect),length(delayVect));
zeromat = zeros(length(incVect),length(delayVect));
prealloc = struct('finalsize',zeromat,'peaktime',zeromat,'peaksize',zeromat,'finaltime',zeromat,'growth',zeromat,'avsused',zeromat,'avsatpeak',zeromat);
dynalloc = struct('finalsize',zeromat,'peaktime',zeromat,'peaksize',zeromat,'finaltime',zeromat,'growth',zeromat,'avsused',zeromat,'avsatpeak',zeromat);

stateList=[];
houseSizes=[];
recNum=[];
for i=1:length(var.pi_k)
	[Q, t]=genQ(var,i,0);
	stateList=[stateList, t];
	houseSizes=[houseSizes, i*ones(1,length(t))];
	recNum=[recNum, i-(t(1,:)+t(2,:)+t(3,:))];
end
pop=sum(stateList(1:3,:))+recNum;
for i=1:length(incVect)
	for j=1:length(delayVect)
		fprintf('Inc: %d Delay: %.4f\n',incVect(i),delayVect(j));
		fprintf('Calculating pre-allocation: \n');
		%Pre-allocation: (independent of delay time)
		%var.maxAV=avVect(i);
		
		var.initialAV=incVect(i);
		var.kappa=delayVect(j);
		%var.maxAV=100000;
		vart=var;
		vart.maxAV=var.initialAV;
		
        var.maxAV=incVect(i);
		var.phi_k=calcProportions(var);
        
		%var.initialAV=incVect(i);
		fprintf('Initial Conditions: ');
		[ic, prealloc.growth(i,j)]=init_conds_mh(var);
		fprintf('Calculating DE''s: ');
		[tout, yout]=selfConstMH(var,ic,[0 100]);
		%[tout, yout, prealloc.growth(i,j)]=solveDE(var,incVect(i),[0 100]);
		
		fprintf('Done! \nDoing statistics... ');
		
		%Statistics:
		prealloc.finalsize(i,j)=recNum*(var.N.*yout(end,:))';
		
		infAtT=var.N*stateList(3,:)*(yout');
		[prealloc.peaksize(i,j), peakIndex]=max(infAtT);
		prealloc.peaktime(i,j)=tout(peakIndex);
		
		prealloc.finaltime(i,j)=tout(end);
		
		prealloc.avsused(i,j)=pop*(var.N.*((stateList(4,:)~=0).*yout(end,:)))';
		prealloc.avsatpeak(i,j)=pop*(var.N.*((stateList(4,:)~=0).*yout(peakIndex,:)))';
		preallocavs{i,j}=(stateList(4,:)~=0)*yout'*var.N*3;
		pret{i,j}=tout;
		fprintf('Done! \n');
		
		
		%Dynamic Allocation:
		fprintf('Calculating dyn-allocation: \n');
		var.kappa=delayVect(j);
		var.phi_k=zeros(1,length(var.pi_k));
		
		%fprintf('Initial Conditions: ');
		[ic, dynalloc.growth(i,j)]=init_conds_mh(var);
		fprintf('Calculating DE''s: ');
		[tout, yout]=selfConstMH(var,ic,[0 100]);
		%[tout, yout, dynalloc.growth(i,j)]=solveDE(var,incVect(i),[0 100]);
		
		%Statistics:
		fprintf('Done! \nDoing statistics...');
		dynalloc.finalsize(i,j)=recNum*(var.N.*yout(end,:))';
		
		infAtT=var.N*stateList(3,:)*(yout');
		[dynalloc.peaksize(i,j), peakIndex]=max(infAtT);
		dynalloc.peaktime(i,j)=tout(peakIndex);
		
		dynalloc.finaltime(i,j)=tout(end);
		
		dynalloc.avsused(i,j)=pop*(var.N.*((stateList(4,:)~=0).*yout(end,:)))';
		dynalloc.avsatpeak(i,j)=pop*(var.N.*((stateList(4,:)~=0).*yout(peakIndex,:)))';
		dynallocavs{i,j}=(stateList(4,:)~=0)*yout'*var.N*3;
		dynt{i,j}=tout;
		fprintf('Done!\n\n');

		
		
		
	end
end