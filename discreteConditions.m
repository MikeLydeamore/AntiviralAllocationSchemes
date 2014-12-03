function dcond = discreteConditions(var,icond,stateList)
%Determines a discrete set of initial conditions, appropriate for use of
%simulation.


I = var.N*icond*stateList(3,:)';
Idcond=0;
Ndcond=0;
ncond = var.N*icond;
l=length(icond);

dcond = floor(ncond); %Take the whole integer part.
econd = zeros(size(dcond));
rounded = ncond-dcond;

while abs(Idcond-I)>0.5 || Ndcond ~= var.N
    rands = rand(1,l);
	 econd=zeros(size(dcond));
    for i=1:l
        if rands(i) < rounded(i)
            econd(i)=econd(i)+1;
        end
    end
    Idcond = (dcond+econd)*stateList(3,:)';
    Ndcond = sum(dcond+econd);

end

dcond=dcond+econd;
end
