function [tout, yout, gr]=solveDE(var,increment,tspan)

%Generate AV available vector.
avAvailable=var.initialAV:increment:var.maxAV;

time=tspan(1):0.4:tspan(2);

%Pad AV Vector
avAvailable=[avAvailable, var.maxAV*ones(1,length(time)-length(avAvailable))];

var.maxAV=avAvailable(1);
if var.zeta==Inf
	var.phi_k=calcProportions(var);
end
[ic, gr]=init_conds_mh(var);

tout=[];
yout=[];

stateList=[];
houseSizes=[];
recNum=[];
for i=1:length(var.pi_k)
    [Q, t]=genQ(var,i,0);
    stateList=[stateList, t];
    houseSizes=[houseSizes, i*ones(1,length(t))];
    recNum=[recNum, i-(t(1,:)+t(2,:)+t(3,:))];
end

for i=1:length(time)-1
    var.maxAV=avAvailable(i);
    [ttemp, ytemp]=selfConstMH(var,ic,[time(i), time(i+1)]);
    ic=ytemp(end,:);
    if var.zeta==Inf
        %Shift mass around
        ic=shiftMass(var,ic,avAvailable(i+1)-avAvailable(i),houseSizes,stateList);
        %disp(avAvailable(i+1)-avAvailable(i));
    end
    tout=[tout; ttemp];
    yout=[yout; ytemp];
    if (houseSizes.*(stateList(4,:)~=0))*yout(end,:)'*var.N > var.maxAV
        x=2;
    end
    if avAvailable(i+1)==avAvailable(i)
        var.maxAV=avAvailable(i+1);
        if i+1~=length(time)
            [ttemp, ytemp]=selfConstMH(var,ic,[time(i+1) time(end)]);
            tout=[tout; ttemp];
            yout=[yout; ytemp];
        end
        break;
    end
end
end

function out=shiftMass(var,in,amountToShift,houseSizes,stateList)
out=in;
smallestHouse=8;
%Work out how many AVs are left to give out
while amountToShift>smallestHouse
    maxToGiveOut=amountToShift;
    
    pullOut=find(stateList(4,:)==0);
    v=out;
    v(stateList(4,:)~=0)=0; %Zero out the elements we're not interested in.
	 %v((stateList(1,:)+stateList(2,:)+stateList(3,:))==houseSizes)=0; %New - added 19/6 - ensures we don't give antivirals to an already completed household.
    u=rand;
    v=cumsum(v./sum(v));
    
    index=find(u<v,1);
    
    houseSize=houseSizes(index);
    while houseSize > amountToShift
        u=rand;
        index=find(u<v,1);
        houseSize=houseSizes(index);
    end
    
    s=stateList(1,index);
    e=stateList(2,index);
    i=stateList(3,index);
    
	 if ((s+e) == houseSizes(index))
		 moveTo=find(s==stateList(1,:)&e==stateList(2,:)&i==stateList(3,:)&4==stateList(4,:));
	 else
		 moveTo=find(s==stateList(1,:)&e==stateList(2,:)&i==stateList(3,:)&1==stateList(4,:));
	 end
    for ii=1:length(moveTo)
        if houseSizes(moveTo(ii))==houseSize
            moveTo=moveTo(ii);
            break;
        end
    end
            

   if stateList(4,index)~=0
       error('Oh dear...');
   end
    out(moveTo)=out(moveTo)+1/var.N;
    out(index)=out(index)-1/var.N;
    amountToShift=amountToShift-houseSize;
    
    temp=out;
    temp(stateList(4,:)~=0)=0;
    smallestHouse=find(temp>eps,1);
    smallestHouse=houseSizes(smallestHouse);
end


end