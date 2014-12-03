%function r_comparison_sim(beta,gamma,sigma,alpha,onDist,zeta,offDist,kappa,rho,tau,eta,pi_k,phi_k,N,iter)
%variables;

testIterations=0;
j=zeros(1,iter);

%tic;
%rs=rStar(var);

%p=toc;
%tic;
%phi_k, noAV]=calcProportions(N,pi_k,maxAV,dosage);
%r=malthusian(var);
%o=toc;
stateList=[];
houseSizes=[];
for k=1:length(var.pi_k)
    [Q, s]=genQ(var,k,0);
    stateList=[stateList, s];
    houseSizes=[houseSizes, k*ones(1,length(s))];
end
%[Q, stateList]=genQHalf(var,k);

%fprintf('r: %.4f - Completed in %.4fs \n',r,o);
%fprintf('R*: %.4f - Completed in %.4fs \n',rs,p);

x=0:0.001:90;
avgInf=zeros(1,length(x));
avgAV=0;
infected=cell(1,iter);
eventTime=cell(1,iter);
noAV=cell(1,iter);
totalInfected=cell(1,iter);
popSize=cell(1,iter);
tic;
%initialConditions=[initialCondition, zeros(1,2*length(initialCondition))];
initialConditions=initialCondition;
parfor i=1:iter
	%testIterations=testIterations+1;
	%fprintf('Iteration: %2d \n',i);
	%tic;
	[infected{i}, eventTime{i}, noAV{i}, totalInfected{i}, popSize{i}]=simulation(var,initialConditions,stateList,houseSizes,var.maxAV);
	%totalTime=eventTime;
	%endtime(i)=totalTime(end);
	%%Work out average
	%avgInf(1)=avgInf(1)+sum(stateList(3,:).*(N*initialConditions)); %start with 20 infected
	%for l=1:length(x)
		
	%	interval=find((x(l)>=totalTime));
%		interval=interval(end);
		%avgInf(l)=avgInf(l)+infected(interval);
	%end
	%avgAV=avgAV+noAV;
	
	%j(i)=toc;
	%fprintf('Number of AV''s allocated: %d | Epidemic end time: %.4f \n Completed in: %.4fs \n',noAV,endtime(i),j(i));
end
toc;
%plot(x,log(avgInf/iter));
%axis([0 50 0 8.5]);
%hold on;
%plot(x,r*x+log(avgInf(2)/iter),'r');
%xlabel('Time');
%ylabel('log(infected)');
%legend('Simulation','r');
%title('log(I) = Time*r');
%fprintf('Average number of infected: %d\n',sum(totalInfected)/iter);
%fprintf('Average number of allocated antivirals: %d\n',avgAV/iter);
