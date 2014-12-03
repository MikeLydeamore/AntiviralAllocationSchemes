varPre=var;
varPre.zeta=Inf;
varPre.phi_k=calcProportions(varPre);
varDyn=var;
varDyn.zeta=1;
popsize=(var.N*var.pi_k)*(1:length(pi_k))';
avVect=1000:1000:popSize;
for i=1:length(avVect)
    fprintf('Currently searching %d\n',avVect(i));
    dynsize=calcFinalSize(varDyn,avVect(i),0);
    p(i)=findprob(varPre,dynsize,avVect(i),1);
end

%plot(avVect,p);