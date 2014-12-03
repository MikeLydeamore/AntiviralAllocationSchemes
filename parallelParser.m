%load bigsim_01.mat
avgAV=0;
x=0:0.001:90;
avgInf=zeros(1,length(x));
for i=1:length(eventTime)
    fprintf('Parsing %d \n',i);
    avgInf(1)=avgInf(1)+infected{i}(1); %start with 20 infected
    for l=1:length(x)
        interval=find((x(l)>=eventTime{i}));
        interval=interval(end);
        avgInf(l)=avgInf(l)+infected{i}(interval);
    end
    avgAV=avgAV+noAV{i};
end