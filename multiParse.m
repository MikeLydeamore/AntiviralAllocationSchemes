x=0:0.01:90;
avgInf=zeros(1,length(x));
figure(1);
hold on;
counter=0;
for i=2:5
	load(['psitest' int2str(i) '.mat'],'eventTime','infected');
	clc;
	fprintf('Onto dataset %d\n',i);
	for j=1:length(eventTime)
		fprintf('Iteration: %d of %d\n',j,length(eventTime));
		avgInf(1)=avgInf(1)+infected{j}(1); %start with 20 infected
		for l=1:length(x)
			interval=find((x(l)>=eventTime{j}));
			interval=interval(end);
			avgInf(l)=avgInf(l)+infected{j}(interval);
		end
		plot(eventTime{j},infected{j});
		counter=counter+1;
	end
end
