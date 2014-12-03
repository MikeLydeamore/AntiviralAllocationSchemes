function rs=rs_prealloc(var)
%Calculates R* for an exp-exp type problem.

rs=0;
var.zeta=1/var.zeta;
var.kappa=1/var.kappa;
for k=1:length(var.pi_k)
	[Q, stateList]=genQ(var,k,0);
	transient=find(diag(Q)~=0);
	
	for i=1:length(transient)
		for j=1:length(transient)
			A1(i,j)=Q(transient(i),transient(j));
		end
		%The external infections happen at rate alpha when there is no antivrials
		%and rate (1-tau)*alpha when there is antivirals.
		if stateList(4,transient(i))==1
			b1(i)=(1-var.tau)*var.alpha*stateList(3,transient(i));
		else
			b1(i)=var.alpha*stateList(3,transient(i));
		end
	end
	
	x=A1\-b1'; %Ax+b = 0
	
	%Want state (k-1, 1, 0, 4)
	initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0&stateList(4,:)==4);
	
	term1=x'*(initialState==transient);
	
	%No preallocation sections
	[Q, stateList]=genQHalf(var,k);
	transient=find(diag(Q)~=0);
	
	for i=1:length(transient)
		for j=1:length(transient)
			A2(i,j)=Q(transient(i),transient(j));
		end
		%The external infections happen at rate alpha
			b2(i)=var.alpha*stateList(3,transient(i));
	end
	
	x2=A2\-b2'; %Ax+b = 0
	
	%Want state (k-1, 1, 0)
	initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==1&stateList(3,:)==0);
	term2=x2'*(initialState==transient);
	
	rs=rs+(var.pi_k(k)*( (var.phi_k(k)*term1) + ((1-var.phi_k(k))*term2)) );
end
end