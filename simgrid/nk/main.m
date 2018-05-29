% nkpfsim.m
% Main file for computing the policy functions in small NK model without ZLB
% solved with simulation-based grid
% Written by Takeki Sunakawa
% May 4 2018  @ Guesthouse Universitat Mannheim
% May 29 2018 @ Starbucks Mannheim

clear all;

%% metaparameters
pfmethod = 1;
ngh = 3; % number of gh nodes
tol = 1e-6; % tolerance for main loop
damp = 0.7; % dampling parameter (for old policy functions)
simT = 10000; % length of simulation for Euler error

%% setup for the EDS grid
% number of EDS grid
nv = 25;
% number of simulation
per = 10000;
% degree of global polynomial
polyd = 2;
% interval for sampling
mu = 5;
% threshold of density
crit = 1e-3;

%% parameter values
% these parameters are chosen by me
invnu = 6; % 1/nu = elasticity
gyss = 0.2;
gbar = 1.0/(1.0-gyss);
% from Herbst and Schorfheide p. 98
rrbar = 0.42;
% pibar = 1.0+3.30/400;
% gamma = 1.0+0.52/100;
% beta  = 1.0/(1.0+rrbar/400);
pibar = exp(3.30/400);
gamma = exp(0.52/100);
beta  = 1.0/exp(rrbar/400);

tau = 2.83;
kappa = 0.78;
phi = tau*(invnu-1)/pibar^2/kappa;
psi1 = 1.80;
psi2 = 0.63;
rhor = 0.77;
rhog = 0.98;
rhoz = 0.88;
sigmar = 0.22/100;
sigmag = 0.71/100;
sigmaz = 0.31/100;

disp(' ');
disp(' number of grid points nv=25 and order of polynomial np=2');
% nkpf2(0,25,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
nkpf2(1,25,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
% 
% disp(' ');
% disp(' number of grid points nv=50 and order of polynomial np=2');
% nkpf2(0,50,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
% nkpf2(1,50,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
% 
% disp(' ');
% disp(' number of grid points nv=100 and order of polynomial np=2');
% nkpf2(0,100,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
% nkpf2(1,100,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);


function nkpf2(pfmethod,nv,per,polyd,mu,crit,ngh,tol,damp,simT,...
    pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz)

tic;

%% setups
peaflag = 1; % =1; expectation term as a whole is approximated (only with pfmethod=0)

% number of variables
ns = 4;
ne = ns-1;
nghe = ngh^ne;

% steady state values
css = (1-1/invnu)^(1/tau);
piss = pibar;
rnss = gamma*pibar/beta;
yss = gbar*css;
fcss = beta*css^(-tau)/gamma/piss;
fpss = beta*phi*css^(-tau)*yss*(piss-pibar)*piss;

% Gaussian-Hermite quadrature
% xz is abscissa and wz is weight for eps
[xg,wg] = qnwnorm(ngh,0,sigmag);
[xz,wz] = qnwnorm(ngh,0,sigmaz);
[xr,wr] = qnwnorm(ngh,0,sigmar);

wmat = zeros(nghe,3);
xmat = zeros(nghe,3);
for ir = 1:ngh

    for iz = 1:ngh

        for ig = 1:ngh

            index = ngh^2*(ir-1)+ngh*(iz-1)+ig;
            xmat(index,:) = [xg(ig) xz(iz) xr(ir)];
            wmat(index,:) = [wg(ig) wz(iz) wr(ir)];

        end

    end

end

% initial grid points: from log-linearized solution
rng(1);
shocks = randn(per,ne);
disp('Simulating the ergodic distribution');
[sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = ...
    calcsim_init(shocks,invnu,gbar,pibar,gamma,beta,tau,kappa,rhor,psi1,psi2,rhog,rhoz,sigmar,sigmag,sigmaz);

disp('Generating the EDS grid');
sim(1:per,1) = sim_rn;
sim(1:per,2) = sim_g;
sim(1:per,3) = sim_z;
sim(1:per,4) = sim_r;
[sf eps] = eds_grid(sim,mu,crit,nv);
nv = size(sf,1);

h = figure;
subplot(231);
plot(sim_rn,sim_g,'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,2),'kx','MarkerSize',7.5);
hold off;
subplot(232);
plot(sim_rn,sim_z,'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,3),'kx','MarkerSize',7.5);
hold off;
subplot(233);
plot(sim_rn,sim_r,'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,4),'kx','MarkerSize',7.5);
hold off;
drawnow;

%% outer loop
conv = 0;
dist = 1d+4;
iterout = 0;

%while (conv==0) 
for iterout = 1:2
    
%    iterout = iterout + 1;
    % convergence criteria: See JMM
    if (dist<eps*2); conv = 1; end;

    % set up the basis functions
    X = makebas4(sf,polyd);
    % LS-SVD
    [U,S,V] = svd(X,0);
    bbtinv = V*inv(S)*U';
    
    % initial values
    cvec0  = css*ones(nv,1);
    cvec1  = zeros(nv,1);
    pivec0 = piss*ones(nv,1);
    pivec1 = zeros(nv,1);
    rnvec0 = rnss*ones(nv,1);
    rnvec1 = zeros(nv,1);
    yvec0  = yss*ones(nv,1);
    yvec1  = zeros(nv,1);

    fcvec0 = fcss*ones(nv,1);
    fcvec1 = zeros(nv,1);
    fpvec0 = fpss*ones(nv,1);
    fpvec1 = zeros(nv,1);

    %% inner loop
    disp('Solving for the policy functions');
    diff = 1e+4;
    iter = 0;

    while ((diff>tol) && (iter<1000))

        % fitting polynomials
        coefc  = bbtinv*cvec0;
        coefpi = bbtinv*pivec0;
        coeffc = bbtinv*fcvec0;
        coeffp = bbtinv*fpvec0;

        for i=1:nv

            rnpast = sf(i,1);
            gnow = sf(i,2);
            znow = sf(i,3);
            rnow = sf(i,4);
            ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);

            if (pfmethod==0)

                % solve nonlinear equations for c and pi
                x0 = [cvec0(i) pivec0(i)]';
                if (peaflag==1)
                   [c0 pi0] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,polyd,xmat,wmat,tol,...
                       pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                else
                   [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,tol,...
                       pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                end
                y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
                rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);            

                % expectation term (used when peaflag==1)
                fc0 = beta*c0^(-tau)/(gamma*exp(znow))/pi0;
                fp0 = beta*phi*c0^(-tau)*y0*(pi0-pibar)*pi0;
                
            elseif (pfmethod==1)
                
                fc0 = fcvec0(i);
                fp0 = fpvec0(i);
                [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi);

                y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
                rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);

                % update the expectation term f with interpolation
                fc0 = 0.0;
                fp0 = 0.0;
                for ighe=1:nghe

                    % next period's fc and fp (obtained by interpolation)                
                    gp = rhog*gnow + xmat(ighe,1);
                    zp = rhoz*znow + xmat(ighe,2);
                    rp = xmat(ighe,3);
                    fcp = makebas4([rn0 gp zp rp],polyd)*coeffc;
                    fpp = makebas4([rn0 gp zp rp],polyd)*coeffp;

                    % next period's c and pi (obtained by next period's fc and fp)                                
                    [cp pip] = pf1(fcp,fpp,pibar,invnu,tau,phi);
                    yp = c0/(1/gbar/exp(gp) - phi/2*(pip-pibar)^2);

                    % current period's fc and fp                
                    fcx = cp^(-tau)*rn0/(gamma*exp(zp))/pip;
                    fpx = cp^(-tau)*yp*(pip-pibar)*pip/y0;
                    weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
                    fc0 = fc0 + weight*beta*fcx;
                    fp0 = fp0 + weight*beta*phi*fpx;

                end

            end

            cvec1(i)  = c0;
            pivec1(i) = pi0;
            yvec1(i)  = y0;
            rnvec1(i) = rn0;
            fcvec1(i) = fc0;
            fpvec1(i) = fp0;

        end

        % calculate the norm between the old and new policy functions
        diffc = max(abs(cvec1-cvec0));
        diffp = max(abs(pivec1-pivec0));
        diffr = max(abs(rnvec1-rnvec0));
        diffy = max(abs(yvec1-yvec0));
        diff  = max([diffc diffp diffr diffy]);

        % update the policy functions
        cvec0  = damp*cvec0 + (1.0-damp)*cvec1;
        pivec0 = damp*pivec0 + (1.0-damp)*pivec1;
        rnvec0 = damp*rnvec0 + (1.0-damp)*rnvec1;
        yvec0  = damp*yvec0 + (1.0-damp)*yvec1;
        fcvec0 = damp*fcvec0 + (1.0-damp)*fcvec1;
        fpvec0 = damp*fpvec0 + (1.0-damp)*fpvec1;

        % counter for inner-loop iterations
        iter = iter + 1;

%        disp([iter diffc diffp diffy]);

    end

    % generate new EDS grid
    disp('Simulating the ergodic distribution');
    [sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = ...
        calcsim(shocks,coefc,coefpi,coeffc,coeffp,pfmethod,peaflag,polyd,xmat,wmat,tol,...
        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss,sigmar,sigmag,sigmaz);

    disp('Generating the EDS grid');
    sim(1:per,1) = sim_rn;
    sim(1:per,2) = sim_g;
    sim(1:per,3) = sim_z;
    sim(1:per,4) = sim_r;
    [sf_up eps] = eds_grid(sim,mu,crit,nv);
    nv_up = size(sf_up,1);
    
    % calculate the norm between the old and new EDS grid
    dist = 0;
    for iv = 1:nv_up

        dist_up = min(sum((ones(nv,1)*sf_up(iv,:)-sf).^2,2).^.5);
        dist = max(dist,dist_up);

    end

    disp(' ');
    disp(sprintf('AT ITERATION ON GRIDS  = %5d', iterout));
    disp('----------------------------------------');
    disp(sprintf('Distance between new and old EDS grids: %1.5f',dist));
    disp(sprintf('Criteria of convergence from EDS grids: %1.5f',eps*2));
    disp('----------------------------------------');

    % update EDS grid
    sf = sf_up;
    nv = nv_up;

    figure(h);
    subplot(231);
    plot(sim_g,sim_rn,'k.','Color',[.5 .5 .5]);
    hold on;
    plot(sf(:,2),sf(:,1),'kx','MarkerSize',7.5);
    ylabel('R_{-1}^*'); xlabel('g')
    hold off;
    subplot(232);
    plot(sim_z,sim_rn,'k.','Color',[.5 .5 .5]);
    hold on;
    plot(sf(:,3),sf(:,1),'kx','MarkerSize',7.5);
    xlabel('z')
    hold off;
    subplot(233);
    plot(sim_r,sim_rn,'k.','Color',[.5 .5 .5]);
    hold on;
    plot(sf(:,4),sf(:,1),'kx','MarkerSize',7.5);
    xlabel('\epsilon_R')
    hold off;
    drawnow;

end

t = toc;

%% Euler errors
drop = floor(0.05*simT);
simTT = simT + drop;

rnvec = zeros(simTT,1);
yvec = zeros(simTT,1);
cvec = zeros(simTT,1);
pivec = zeros(simTT,1);
gvec = zeros(simTT,1);
zvec = zeros(simTT,1);
rvec = zeros(simTT,1);
rnvec(1) = rnss;
rng(0);

for time = 1:simTT-1
    
    rnpast = rnvec(time);
    gnow = gvec(time);
    znow = zvec(time);
    rnow = rvec(time);
    ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);

    if (pfmethod==0)

        c0 = makebas4([rnpast gnow znow rnow],polyd)*coefc;
        pi0 = makebas4([rnpast gnow znow rnow],polyd)*coefpi;
        
    elseif (pfmethod==1)

        fc0 = makebas4([rnpast gnow znow rnow],polyd)*coeffc;
        fp0 = makebas4([rnpast gnow znow rnow],polyd)*coeffp;        
        [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi);
        
    end

    y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
    rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);
    
    % euler errors
    fc0 = 0.0;
    fp0 = 0.0;
    for ighe=1:nghe

        gp = rhog*gnow + xmat(ighe,1);
        zp = rhoz*znow + xmat(ighe,2);
        rp = xmat(ighe,3);

        if (pfmethod==0)

            cp = makebas4([rn0 gp zp rp],polyd)*coefc;
            pip = makebas4([rn0 gp zp rp],polyd)*coefpi;
            
        elseif (pfmethod==1)

            fcp = makebas4([rn0 gp zp rp],polyd)*coeffc;
            fpp = makebas4([rn0 gp zp rp],polyd)*coeffp;
            [cp pip] = pf1(fcp,fpp,pibar,invnu,tau,phi);
            
        end
        
        yp = cp/(1/gbar/exp(gp) - phi/2*(pip-pibar)^2);
        fcx = beta*cp^(-tau)*rn0/(gamma*exp(zp))/pip;
        fpx = beta*phi*cp^(-tau)*yp*(pip-pibar)*pip;
        
        weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
        fc0 = fc0 + weight*fcx;
        fp0 = fp0 + weight*fpx;

    end
    
    evec(time+1,1) = 1 - (c0^tau)*fc0;
    LHS2 = ( (1-invnu)+invnu*c0^tau - phi*(pi0-pibar)*(pi0-.5*invnu*(pi0-pibar)) )*c0^(-tau)*y0;
    evec(time+1,2) = LHS2 + fp0;
    
    rnvec(time+1) = rn0;
    yvec(time+1)  = y0;
    cvec(time+1)  = c0;
    pivec(time+1) = pi0;
    gvec(time+1)  = rhog*gnow + sigmag*randn;
    zvec(time+1)  = rhoz*znow + sigmaz*randn;
    rvec(time+1)  = sigmar*randn;
    
%    disp([c0 pi0 time]);
    
end

rnvec = rnvec(drop+1:simTT);
dyvec = yvec(drop+1:simTT)./yvec(drop:simTT-1);
cvec  = cvec(drop+1:simTT);
pivec = pivec(drop+1:simTT);

disp([std(log(dyvec)*100) std(log(pivec)*400) std(log(rnvec)*400)]);
disp([log10(mean(abs(evec))) log10(max(abs(evec))) t iter]);
%disp([log10(mean(abs(evec(:,1)))) log10(max(abs(evec(:,1)))) t]);

end


function [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2PV',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

end


function [c0 pi0] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,polyd,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2PVf',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,polyd,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

end


function [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi)

c0 = fc0^(-1/tau);
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end