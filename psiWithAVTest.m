avVect = 1000:1000:100000;
prop = zeros(1,length(avVect));

for i=1:length(avVect)
	var.maxAV=avVect(i);
	prop(i) = psiZero(var);
end
	