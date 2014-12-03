function av=findAV(var)

av=fzero(@(a) CFS(var,a),65000);


end

function fs = CFS(var,av)
	var.maxAV=av;
	varpre=var;

	varpre.zeta=Inf;
	varpre.phi_k=calcProportions(varpre);
	ic1=init_conds_mh(varpre);
	
	var.zeta=0.5;
	var.phi_k=zeros(1,length(var.pi_k));
	ic2=init_conds_mh(var);
	
	fprintf('Testing %d\n',av);
	
	fs = calcFinalSize(varpre,av,1,ic1)-calcFinalSize(var,av,0,ic2);
	fprintf('Difference: %f\n',fs);

end