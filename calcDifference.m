function r = calcDifference(b,var)

	var.beta=b;
    fprintf('Testing Beta: %.4f\n',b);
	
	r = calcFinalSize(var,0,1,init_conds_mh(var));


end