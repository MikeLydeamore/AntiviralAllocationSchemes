function s=malthusian_f(r,var)

%Zeta is the expected delay to antiviral arrival, meaning that the delay is
%exp~(1/zeta).
%Kappa is the expected duration of antivirals.

%Break up into 3 Integrals
s=0;
var.zeta = 1/var.zeta;
for k=1:length(var.pi_k)
    %Integral 1
    %0-T where T is ~exp(1/zeta).
    %Use the Q matrix with absorption after antivirals introduced.
    [Q, stateList]=genQ(var,k,1);
    transient=find(diag(Q)~=0);
    
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    for i=1:length(transient)
        for j=1:length(transient)
            A1(i,j)=Q(transient(i),transient(j));
        end
        %No antivirals here.
        b1(i)=var.alpha*stateList(3,transient(i));
    end
    
    x=A1\-b1'; %Ax+b = 0
    
    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1&stateList(4,:)==0);
    integral1=x'*(initialState==transient);
    
    %Integral 2
    
    %First term: p(d)zexp(-d(z+r)) dd
    noAVvar=var;
    noAVvar.tau=0;
    noAVvar.rho=0;
    noAVvar.eta=0;
    [Q, stateList]=genQHalf(noAVvar,k);
    p0=zeros(1,length(stateList));
    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1);
    p0(initialState)=1;
    first_term=p0*(var.zeta)*inv(-(Q-(var.zeta+r)*eye(size(Q))));
    
    
    %Second term - Integral from 0-kappa, using phiv with the reduced Q, as
    %normal.
    [Q, stateList]=genQHalf(var,k);
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    for i=1:length(transient)
        for j=1:length(transient)
            A2(i,j)=Q(transient(i),transient(j));
        end
        %There is always antivirals for this integral.
        b2(i)=(1-var.tau)*var.alpha*stateList(3,transient(i));
    end
    
    psi0=zeros(length(b2),1);
    second_term = phiv(var.kappa,A2,b2',psi0);
    
    integral2=first_term(transient)*second_term;
    
    %Integral 3
    %First term:
    [Q1, stateList]=genQHalf(noAVvar,k);
    [Q2, stateList]=genQHalf(var,k);
    
    p0=zeros(1,length(stateList));
    initialState=find(stateList(1,:)==(k-1)&stateList(2,:)==0&stateList(3,:)==1);
    p0(initialState)=1;
    
%    first_term=expv(-var.kappa,(r*eye(size(Q2))-Q2)',p0*(var.zeta)*inv((var.zeta+r)*eye(size(Q1))-Q1));
	first_term = var.zeta*p0*expm(-var.kappa*(r*eye(size(Q2))-Q2))*inv((var.zeta+r)*eye(size(Q1))-Q1);
    
    %Second term - Integral from 0-Inf of unrreduced system.
    [Q, stateList]=genQHalf(noAVvar,k);
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    for i=1:length(transient)
        for j=1:length(transient)
            A3(i,j)=Q(transient(i),transient(j));
        end
        %There is never antivirals for this integral.
        b3(i)=var.alpha*stateList(3,transient(i));
    end
    
    second_term=A3\-b3';
    
    integral3=first_term(transient)*second_term;
    
    noPreAllocation=integral1+integral2+integral3;
    
    %Incorporate pre-allocated households
    
    %First integral - 0-Kappa (1-tau)
    %p(0) = (k-1,0,1)
    [Q, stateList]=genQHalf(var,k);
    
    p0=zeros(1,length(stateList));
    p0(stateList(1,:)==k-1&stateList(2,:)==0&stateList(3,:)==1)=1;
    
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    for i=1:length(transient)
        for j=1:length(transient)
            A4(i,j)=Q(transient(i),transient(j));
        end
        %There is always antivirals for this integral.
        b4(i)=(1-var.tau)*var.alpha*stateList(3,transient(i));
    end
    
    psi0=zeros(length(b2),1);
    integral = phiv(var.kappa,A4,b4',psi0);
    
    integral1=p0(transient)*integral;
    
    %Second integral - Kappa-Inf, unreduced.
    %Get p(Kappa) = p(0)exp(-Q*Kappa) with reduced Q.
    [Q, stateList]=genQHalf(var,k);
    p0=zeros(1,length(stateList));
    p0(stateList(1,:)==k-1&stateList(2,:)==0&stateList(3,:)==1)=1;
    pK=mexpv(var.kappa,Q',p0);
    
    %Now use the unreduced Q for the remainder
    [Q, stateList]=genQHalf(noAVvar,k);
    transient=find(diag(Q)~=0);
    for i=1:length(transient)
        Q(transient(i),transient(i))=Q(transient(i),transient(i))-r; %Add discounting rate, -(qii + r) = -qii - r
    end
    for i=1:length(transient)
        for j=1:length(transient)
            A5(i,j)=Q(transient(i),transient(j));
        end
        %There is always antivirals for this integral.
        b5(i)=var.alpha*stateList(3,transient(i));
    end
    
    integral=A5\-b5';
    
    integral2=pK(transient)'*integral*exp(-r*var.kappa);
    
    preAllocation=integral1+integral2;
    
    s=s+(var.pi_k(k)*((1-var.phi_k(k))*noPreAllocation+var.phi_k(k)*preAllocation));
end
s=s-((r+var.sigma)/var.sigma);