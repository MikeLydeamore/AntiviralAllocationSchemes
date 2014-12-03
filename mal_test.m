function s = mal_test(r,var)

s=0;
var.zeta = 1/var.zeta;
var.kappa = 1/var.kappa;
for k=1:length(var.pi_k)
    [Q, stateList]=genQ(var,k,0);
	%[Q, stateList]=genQHalf(var,k);
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    
    for i=1:length(transient)
        for j=1:length(transient)
            A(i,j)=Q(transient(i),transient(j));
        end
        %The external infections happen at rate alpha when there is no antivrials
        %and rate (1-tau)*alpha when there is antivirals.
        if stateList(4,transient(i))==1
            b(i)=(1-var.tau)*var.alpha*stateList(3,transient(i));
        else
            b(i)=var.alpha*stateList(3,transient(i));
        end
    end
    
    x=A\-b'; %Ax+b = 0
    
    %Want state (k-1, 1, 0, 0)
    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0&stateList(4,:)==0);
	% initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0);
    
	initialState2=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0&stateList(4,:)==4);
    s=s+(var.pi_k(k)*((1-var.phi_k(k))*x'*(initialState==transient)+var.phi_k(k)*x'*(initialState2==transient)));
end
s=s-1;
end