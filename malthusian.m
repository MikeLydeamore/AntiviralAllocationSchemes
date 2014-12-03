function r=malthusian(var)
%Solves for the malthusian r for an SEIR system without antivirals
%[Q, stateList]=genQ(beta,gamma,sigma,ro,tau,zeta,pi_k);
%Note that zeta and kappa are the expected delay/durations respectively.


if strcmp(var.onDist,'exp') && strcmp(var.offDist,'exp')
	f=@(r) mal_test(r,var);
	r=fzero(f,[0 5]);
    
elseif strcmp(var.onDist,'exp') && strcmp(var.offDist,'const')
    f=@(r) mal_expconst(r,var);
    r=fzero(f,[0.01 5]);
    
else
	error('Unknown type of delay or active time distribution, or I''ve ruined the code. Probably the latter');
end

end