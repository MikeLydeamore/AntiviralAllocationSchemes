function p = psiZero(var)

var.phi_k=zeros(1,length(var.pi_k));
ic=init_conds_mh(var);
dynsize = calcFinalSize(var,var.maxAV,0,ic);
var.zeta=Inf;
var.phi_k=calcProportions(var);
var.psi=0;
ic=init_conds_mh(var);
zeroit = @(psi) (dynsize-calcSizeWithPsi(var,psi,ic));
if zeroit(0) < 0 %Dynamic already better
	zeroPsi=0;
else
	zeroPsi = fzero(zeroit,[0 1]);
end

%Calc # of houses who take incorrectly:
var.psi=zeroPsi;
ic=init_conds_mh(var);
[tout, yout]=selfConstMH(var,ic,[0 90]);

stateList=[];
houseSizes=[];
for k=1:length(var.pi_k)
    [Q, t]=genQ(var,k,0);
    stateList=[stateList, t];
    houseSizes=[houseSizes,k*ones(1,length(t))];
end
falseStates = find((stateList(1,:)+stateList(2,:))==houseSizes&stateList(4,:)==1);
falseOverT = yout(:,falseStates);
falseOverT = sum(falseOverT,2);
dt = tout(2:end)-tout(1:end-1);
p = falseOverT(1:end-1)'*dt;

end

function fes = calcSizeWithPsi(var,psi,ic)

var.psi=psi;
fes = calcFinalSize(var,var.maxAV,1,ic);
disp(psi);

end