testIterations=0;
j=zeros(1,iter);

stateList=[];
houseSizes=[];
for k=1:length(var.pi_k)
    [Q, s]=genQ(var,k,0);
    stateList=[stateList, s];
    houseSizes=[houseSizes, k*ones(1,length(s))];
end

x=0:0.001:90;
avgInf=zeros(1,length(x));
avgAV=0;
infected=cell(1,iter);
eventTime=cell(1,iter);
noAV=cell(1,iter);
totalInfected=cell(1,iter);
popSize=cell(1,iter);
tic;

initialConditions=initialCondition;
parfor i=1:iter

	[infected{i}, eventTime{i}, noAV{i}, totalInfected{i}, popSize{i}]=simulation(var,initialConditions,stateList,houseSizes,var.maxAV);
	
end
toc;

