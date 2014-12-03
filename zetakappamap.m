counter1=1;
counter2=1;
%variables;

for z=0.01:0.1:5
	var.zeta=z;
	for k=0.01:0.1:5
		fprintf('Zeta: %.2f Kappa: %.2f\n',z,k);
		var.kappa=k;
		heater(counter1,counter2)=malthusian(var);
		counter2=counter2+1;
	end
	counter1=counter1+1;
	counter2=1;
end
	