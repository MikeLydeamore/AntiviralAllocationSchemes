function final = finalsizeSIRode(beta,gamma,N,S0)

rho = gamma*N/beta;

final = N+rho*lambertw(-(S0*exp(-N/rho))/rho);

end