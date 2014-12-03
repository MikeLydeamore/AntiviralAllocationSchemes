avVect = 0:2000:100000;
%avVect=80000;
hybridPercent = 0:5:100;

finalsize = zeros(length(avVect),length(hybridPercent));

stateList=[];
recNum=[];
for k=1:length(var.pi_k)
    [Q, t]=genQ(var,k,0);
    stateList=[stateList, t];
    recNum=[recNum, k-(t(1,:)+t(2,:)+t(3,:))];
end
pop=sum(stateList(1:3,:))+recNum;

for i=1:length(avVect)
    for j=1:length(hybridPercent)
        fprintf('AV:s %d\n   Percentage: %d\n',avVect(i),hybridPercent(j));
        var.maxAV=floor(hybridPercent(j)/100 * avVect(i));
        var.phi_k=calcProportions(var);
        var.maxAV=avVect(i);
        var.zeta=0.5;
        
        ic=init_conds_mh(var);
        [tout, yout]=selfConstMH(var,ic,[0 90]);
        finalsize(i,j)=var.N*recNum*yout(end,:)';
    end
end

dynsize=zeros(length(avVect),length(hybridPercent));
presize=dynsize;

for i=1:length(avVect)
    var.maxAV=avVect(i);
    var.phi_k=zeros(1,length(var.phi_k));
    var.zeta=0.5;
    ic = init_conds_mh(var);
    [tout, yout]=selfConstMH(var,ic,[0 90]);
    dynsize(i,:) = ones(1,length(hybridPercent)).*(var.N*recNum*yout(end,:)');
    
    var.zeta=Inf;
    var.phi_k=calcProportions(var);
    ic=init_conds_mh(var);
    [tout, yout]=selfConstMH(var,ic,[0 90]);
    presize(i,:) = ones(1,length(hybridPercent)).*(var.N*recNum*yout(end,:)');
end