function status = destatus(t,y,flag,var,stateList)
if strcmp(flag,'init')==0 && strcmp(flag,'done')==0
    
    [latestt, index]=max(t);
    inf = var.N*(stateList(3,:)+stateList(2,:))*y(:,index);
    if inf < 5
        status=1;
    else
        status=0;
    end
else
    status=1;
end