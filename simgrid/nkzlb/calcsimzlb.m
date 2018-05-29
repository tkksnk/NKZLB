function [sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = calcsimzlb(shocks,coefcn,coefpin,coefrnn,coeffcn,coeffpn,coefcb,coefpib,coeffcb,coeffpb,pfmethod,peaflag,polyd,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss,sigmar,sigmag,sigmaz)

per = size(shocks,1);

% steady state values
css = (1-1/invnu)^(1/tau);
piss = pibar;
rnss = gamma*pibar/beta;
yss = gbar*css;

sim_c  = css*ones(per,1);
sim_pi = piss*ones(per,1);
sim_rn = rnss*ones(per,1);
sim_y  = yss*ones(per,1);
sim_g  = zeros(per,1);
sim_z  = zeros(per,1);
sim_r  = zeros(per,1);

for i_sim = 1:per-1
    
    rnpast = sim_rn(i_sim);
    rnow = sigmar*shocks(i_sim,1);
    gnow = rhog*sim_g(i_sim) + sigmag*shocks(i_sim,2);
    znow = rhoz*sim_z(i_sim) + sigmaz*shocks(i_sim,3);
    ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);

    if (pfmethod==0)

        rn0 = makebas4([rnpast gnow znow rnow],polyd)*coefrnn;
        if (rn0>=1.0)    
            c0 = makebas4([rnpast gnow znow rnow],polyd)*coefcn;
            pi0 = makebas4([rnpast gnow znow rnow],polyd)*coefpin;
        else
            c0 = makebas4([rnpast gnow znow rnow],polyd)*coefcb;
            pi0 = makebas4([rnpast gnow znow rnow],polyd)*coefpib;
        end
        
    elseif (pfmethod==1)

        rn0 = makebas4([rnpast gnow znow rnow],polyd)*coefrnn;
        if (rn0>=1.0)    
            fc0 = makebas4([rnpast gnow znow rnow],polyd)*coeffcn;
            fp0 = makebas4([rnpast gnow znow rnow],polyd)*coeffpn;
        else
            fc0 = makebas4([rnpast gnow znow rnow],polyd)*coeffcb;
            fp0 = makebas4([rnpast gnow znow rnow],polyd)*coeffpb;
        end        
        
        [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi);
        
    end

    y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
    rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);            

    sim_c(i_sim+1)  = c0;
    sim_pi(i_sim+1) = pi0;
    sim_y(i_sim+1)  = y0;
    sim_rn(i_sim+1) = rn0;
    sim_g(i_sim+1)  = gnow;
    sim_z(i_sim+1)  = znow;
    sim_r(i_sim+1)  = rnow;

end

end


function [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi)

c0 = fc0^(-1/tau);
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end