function [infected, eventTime, noAV, totalInfected, popSize]=simulation(var,initCond,stateList,houseSizes,varargin)
%Simulates the Markov chain where antivirals are introduced.
%beta: Internal Infection Rate
%gamma: Recovery Rate
%sigma: Progression Rate
%alpha: External Infection Rate
%onDist: Distribution of introduction rate (const or exp)
%zeta: Expected delay until antiviral intervention
%offDist: Distribution of removal rate (const or exp)
%kappa: Expected duration of antivirals
%Rho: Reduction in susecptibility due to antivirals
%Tau: Reduction in infectivity due to antivirals
%Eta: Reduction in infectious period due to antivirals
%pi_k: Distribution of household sizes
%phi_k: Distribution of pre-allocation to household sizes (Must be the same length as pi_k)
%N: Number of households
%
%Note that the reduction in infectious period is equivalent to a faster
%recovery rate.


	if isempty(varargin)
		[infected, eventTime, noAV, totalInfected, popSize]=sim_run(beta,gamma,sigma,alpha,zeta,kappa,rho,tau,eta,pi_k,phi_k,N);
	else
		if strcmp(var.onDist,'exp') && strcmp(var.offDist,'exp')
			[infected, eventTime, noAV, totalInfected, popSize]=sim_run_expexp(var,initCond,stateList,houseSizes);
		end
		%[infected, eventTime, noAV, totalInfected, popSize]=sim_run(beta,gamma,sigma,alpha,zeta,kappa,rho,tau,eta,pi_k,N,dosage,varargin);
	end

	
end
