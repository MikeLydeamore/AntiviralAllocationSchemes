function [prop, noAV]=calcProportions(var)

%Calculate household sizes
householdSizes=round(var.h_k(1:end-1).*var.N);
householdSizes(end+1)=var.N-sum(householdSizes);
householdSizesOrig=householdSizes;

%Initialise proportions
prop=zeros(1,length(var.h_k));

availableHouses=1:var.N;
counter=1;
noAV=0;
smallestHouseSize=find(var.h_k>0,1);
prop=zeros(1,length(householdSizes));
while abs(noAV-var.maxAV)>smallestHouseSize && counter<var.N
    
    weighting=cumsum(householdSizes./sum(householdSizes));
    u=rand;
    index=find(u<weighting,1);

    while noAV+var.dosage*index > var.maxAV || householdSizes(index)-1<0
        u=rand;
        index=find(u<weighting,1);
    end
    noAV=noAV+var.dosage*index;
    prop(index)=prop(index)+1/householdSizesOrig(index);
        
    householdSizes(index)=householdSizes(index)-1;
    counter=counter+1;
    if householdSizes(index)==0
        smallestHouseSize=find(householdSizes>0,1);
    end
end

end