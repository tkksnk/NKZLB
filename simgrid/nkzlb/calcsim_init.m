function [sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = calcsim_init(shocks,invnu,gbar,pibar,gamma,beta,tau,kappa,rhor,psi1,psi2,rhog,rhoz,sigmar,sigmag,sigmaz)

n = 6;
nk = 3;
% rn g z c(+1) pi(+1) y(+1)
A = zeros(n,n);
B = zeros(n,n);

% rn - tau*(c(+1)-c) - z(+1) - pi(+1) = 0;
A(1,1) = -1;
A(1,3) = 1;
A(1,4) = tau;
A(1,5) = 1;
B(1,4) = tau;

% -pi + kappa*c + beta*pi(+1) = 0;
A(2,5) = -beta;
B(2,4) = kappa;
B(2,5) = -1;

% -c + y - g = 0;
A(3,2) = 1;
B(3,4) = -1;
B(3,6) = 1;

% -rn + rhor*rn(-1) + (1-rhor)*(psi1*pi+psi2*(y-g)) + eps_r = 0;
A(4,1) = 1;
A(4,2) = (1-rhor)*psi2;
B(4,1) = rhor;
B(4,5) = (1-rhor)*psi1;
B(4,6) = (1-rhor)*psi2;

% -g + rhog*g(-1) + eps_g;
A(5,2) = 1;
B(5,2) = rhog;

% -z + rhoz*z(-1) + eps_z;
A(6,3) = 1;
B(6,3) = rhoz;

[f,p] = solab(A,B,nk);

% POLICY AND TRANSITION FUNCTIONS
%                                       c               pi              rn              y               z               g 
% rn(-1)                     -0.455753       -0.635333        0.440934       -0.455753               0               0
% g(-1)                              0               0               0        0.980000               0        0.980000
% z(-1)                       0.570662        1.023358        0.506359        0.570662        0.880000               0
% eps_z                       0.648480        1.162906        0.575408        0.648480        1.000000               0
% eps_g                              0               0               0        1.000000               0        1.000000
% eps_r                      -0.591887       -0.825108        0.572641       -0.591887               0               0

per = size(shocks,1);
s = zeros(3,3);
s(1,1) = sigmar;
s(2,2) = sigmag;
s(3,3) = sigmaz;

sim_c = zeros(per,1);
sim_pi = zeros(per,1);
sim_rn = zeros(per,1);
sim_y = zeros(per,1);
sim_g = zeros(per,1);
sim_z = zeros(per,1);
sim_r = zeros(per,1);

for i_sim = 1:per-1
    
    rnpast = sim_rn(i_sim);
    gpast = sim_g(i_sim);
    zpast = sim_z(i_sim);
    
    xvec  = [rnpast gpast zpast]';
    yvec  = f*xvec;
    xpvec = p*xvec + s*shocks(i_sim,:)'; %[0 shocks(i_sim,2) 0]';
    
    sim_c(i_sim+1)  = yvec(1);
    sim_pi(i_sim+1) = yvec(2);
    sim_y(i_sim+1)  = yvec(3);
    sim_rn(i_sim+1) = xpvec(1);
    sim_g(i_sim+1)  = xpvec(2);
    sim_z(i_sim+1)  = xpvec(3);
    sim_r(i_sim+1)  = sigmar*shocks(i_sim,1);

end

% steady state values
css = (1-1/invnu)^(1/tau);
piss = pibar;
rnss = gamma*pibar/beta;
yss = gbar*css;

sim_c  = css*exp(sim_c);
sim_pi = piss*exp(sim_pi);
sim_y  = yss*exp(sim_y);
sim_rn = rnss*exp(sim_rn);
% sim_g  = gbar*exp(sim_g);
% sim_z  = gamma*exp(sim_z);