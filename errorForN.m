NVector = [3880 7761 15521 31042];
multiple = [10000 20000 40000 80000];
cmap=['b' 'g' 'm' 'k' 'c' 'r'];
load simwith-3880
divisions = ((1:8)*var.h_k').*NVector;
subplot(1.3,1.3,1);
for supermegacounter=1:length(NVector)

	clear infected tout yout avgInf initialCondition
	%NVect=[3880 7760 15521 23281 31042 38803];
	load(['simwith-' int2str(NVector(supermegacounter))]);
	var.nu=0; var.vaccTime=Inf;	
	[tout, yout]=selfConstMH(var,initialCondition,x);
	
	avgInf = avgInf * multiple(supermegacounter);
	
	plot(tout(1:20:end), abs(((avgInf(1:20:end))-(stateList(3,:)*yout(1:20:end,:)'*var.N))./divisions(supermegacounter)),cmap(supermegacounter));
	drawnow;
	hold on;
    
   % figure(2);
   % plot(tout(1:20:end), abs(((avgInf(1:20:end))-(stateList(3,:)*yout(1:20:end,:)'*var.N)./divisions(supermegacounter))/((avgInf(1:20:end)))),cmap(supermegacounter));
   % drawnow;
   % hold on;
end

xlabel('Time');
ylabel('|Difference|');
legend('10000', '20000','40000','80000');