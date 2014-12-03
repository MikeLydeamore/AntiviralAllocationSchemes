sigma=1;
rho=0.3;
tau=0.3;
onDist='exp';
zeta=0.5;
kappa=1;
offDist='exp';
gamma=1;
eta=0;
k=3;


aus=[1888696 2636361 1240173 1217236 523822 178896 45168 29965];
h_k=aus./sum(aus);
%sudan = [96,145,197,222,244,237,229,195,140,115,80,57,35,21,43]/2056;
%h_k  = sudan;
%pi_k=sudan;
%indo = [5.1,11.1,19.2,23.6,17.9,11.0,5.6,2.9,3.3]/99.7;
%h_k=indo;
%N=sum(pop);

N=round(100000/(h_k*(1:length(h_k))'));
phi_k = zeros(1,length(h_k));

gam=0;
iter=15;
maxAV=71042;
dosage=1;
initialAV=2000;
psi=0;

pi_k = ((1:length(h_k)).*h_k)/(h_k*(1:length(h_k))');

%phi_k=zeros(1,length(pi_k));

var.vaccTime=Inf;
var.nu=0;

%%%Severe
beta = 1.1259;
alpha = 1;

%%%Mild
%beta = 0.9669;
%alpha = 0.8;


%%%% Use this for checking the jacobian
%alpha=1;
%beta=10;
%gamma=100;
%alpha=1000;
var = struct('beta',beta,'gamma',gamma,'sigma',sigma,'alpha',alpha,...
    'eta',eta,'rho',rho,'tau',tau,'psi',psi,'onDist',onDist,'offDist',offDist,'zeta',zeta,...
    'kappa',kappa,'h_k',h_k,'phi_k',phi_k,'pi_k',pi_k,'N',N,'iter',iter,'maxAV',maxAV,...
    'dosage',dosage,'k',k,'initialAV',initialAV);
