NVect=[20000 40000 60000 80000 10000];
tVect = 0:0.001:90;
difference=zeros(5,length(tVect));
for scounter=1:length(NVect)
	clear tout yout infected
	var.N=round(NVect(scounter)/(h_k*(1:length(h_k))'));
	initialCondition = init_conds_mh(var);
	
	r_comparison_sim;
	parallelParser;
	avgInf = (avgInf./var.iter)./NVect(scounter);
	
	[tout, yout]=selfConstMH(var,initialCondition,tVect);
	deInf = stateList(3,:)*yout'*var.N./NVect(scounter);
	clear tout yout infected
	save(['simwith-',num2str(var.N)]);
	
end