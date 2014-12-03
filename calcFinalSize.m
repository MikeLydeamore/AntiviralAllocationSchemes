function fs = calcFinalSize(var,av,preOrNot,ic)
%fprintf('Testing %d\n',av);
stateList=[];
houseSizes=[];
for k=1:length(var.pi_k)
    [Q, t]=genQ(var,k,0);
    stateList=[stateList, t];
    houseSizes=[houseSizes,k*ones(1,length(t))];
end
var.maxAV=av;
if preOrNot==1
%    var.phi_k=calcProportions(var);
%    var.zeta=Inf;
elseif preOrNot==2
    
else
%    var.phi_k=zeros(1,length(var.pi_k));
end

%[ic, gr]=init_conds_mh(var);
[tout, yout]=selfConstMH(var,ic,[0 90]);

rec=houseSizes-(stateList(1,:)+stateList(2,:)+stateList(3,:));
fs=var.N*rec*yout(end,:)';

end