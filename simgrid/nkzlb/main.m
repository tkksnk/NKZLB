% nkpfsimzlb.m
% Main file for computing the policy functions in small NK model with ZLB
% solved with simulation-based grid
% Written by Takeki Sunakawa
% Last updated: May 29 2018

clear all;

%% metaparameters
% pfmethod = 1; % =0: TI, =1: future PEA
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
mu = 5; % 10 with nv=50 does not work
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
nkpf2(0,25,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
% nkpf2(1,25,10000,polyd,mu,crit,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,kappa,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz);
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
maxiterout = 2;
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
disp('Initialization');
rng(1);
shocks = randn(per,ne);
disp(' Simulating the ergodic distribution');
[sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = ...
    calcsim_init(shocks,invnu,gbar,pibar,gamma,beta,tau,kappa,rhor,psi1,psi2,rhog,rhoz,sigmar,sigmag,sigmaz);

% disp(' ');
disp(' Generating the EDS grid');
sim(:,1) = sim_rn;
sim(:,2) = sim_g;
sim(:,3) = sim_z;
sim(:,4) = sim_r;
[sf eps] = eds_grid(sim,mu,crit,nv);
nv = size(sf,1);

h = figure;
subplot(231);
plot(sim(:,1),sim(:,2),'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,2),'kx','MarkerSize',7.5);
hold off;
subplot(232);
plot(sim(:,1),sim(:,3),'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,3),'kx','MarkerSize',7.5);
hold off;
subplot(233);
plot(sim(:,1),sim(:,4),'k.','Color',[.5 .5 .5]);
hold on;
plot(sf(:,1),sf(:,4),'kx','MarkerSize',7.5);
hold off;
drawnow;

%% outer loop
conv = 0;
dist = 1d+4;
iterout = 0;

%while (conv==0) 
for iterout = 1:maxiterout
    
%    iterout = iterout + 1;
    % convergence criteria: See JMM
    if (dist<eps*2); conv = 1; end;

    % set up the basis functions
    X = makebas4(sf,polyd);
    % LS-SVD
    [U,S,V] = svd(X,0);
    bbtinv = V*inv(S)*U';
    
    % initial values
    % NOTE: We use an index-function approach with a pair of policy functions.
    % One assumes the ZLB always binds and the other assumes the ZLB never
    % binds. The next period's policy function is given by a weighted average
    % of the policy function in the ZLB regime and the policy function in the
    % non-ZLB regime with an indicator function. The value of the indicator 
    % function is one when the notional rate is greater than the ZLB, otherwise
    % zero.
    cvec0n  = css*ones(nv,1);
    cvec1n  = zeros(nv,1);
    pivec0n = piss*ones(nv,1);
    pivec1n = zeros(nv,1);
    rnvec0n = rnss*ones(nv,1);
    rnvec1n = zeros(nv,1);
    yvec0n  = yss*ones(nv,1);
    yvec1n  = zeros(nv,1);
    cvec0b  = css*ones(nv,1);
    cvec1b  = zeros(nv,1);
    pivec0b = piss*ones(nv,1);
    pivec1b = zeros(nv,1);
    rnvec0b = rnss*ones(nv,1);
    rnvec1b = zeros(nv,1);
    yvec0b  = yss*ones(nv,1);
    yvec1b  = zeros(nv,1);

    fcvec0n = fcss*ones(nv,1);
    fcvec1n = zeros(nv,1);
    fpvec0n = fpss*ones(nv,1);
    fpvec1n = zeros(nv,1);
    fcvec0b = beta*css^(-tau)/gamma/piss*ones(nv,1);
    fcvec1b = zeros(nv,1);
    fpvec0b = fpss*ones(nv,1);
    fpvec1b = zeros(nv,1);

    %% inner loop
    disp(' ');
    disp(' Solving for the policy functions');
    diff = 1e+4;
    iter = 0;

    while ((diff>tol) && (iter<1000))

        % fitting polynomials
        coefcn  = bbtinv*cvec0n;
        coefpin = bbtinv*pivec0n;
        coefrnn = bbtinv*rnvec0n;
        coeffcn = bbtinv*fcvec0n;
        coeffpn = bbtinv*fpvec0n;
        coefcb  = bbtinv*cvec0b;
        coefpib = bbtinv*pivec0b;
        coeffcb = bbtinv*fcvec0b;
        coeffpb = bbtinv*fpvec0b;

        for i=1:nv

            rnpast = sf(i,1);
            gnow = sf(i,2);
            znow = sf(i,3);
            rnow = sf(i,4);
            ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);

            % time iteration
            if (pfmethod==0)

                % solve nonlinear equations for c and pi
                % in the non-ZLB regime
                x0n = [cvec0n(i) pivec0n(i)]';
                if (peaflag==1)
                    [c0n pi0n] = pf0f(x0n,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,0,polyd,xmat,wmat,tol,...
                        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                else                   
                    [c0n pi0n] = pf0(x0n,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,tol,...
                        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                end
                y0n = c0n/(1/gbar/exp(gnow) - phi/2*(pi0n-pibar)^2);
                rn0n = rnpast^rhor*( rnss*(pi0n/pibar)^psi1*(y0n/ystar)^psi2 )^(1-rhor)*exp(rnow);            

                % in the ZLB regime
                x0b = [cvec0b(i) pivec0b(i)]';

                if (peaflag==1)
                    [c0b pi0b] = pf0f(x0b,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,1,polyd,xmat,wmat,tol,...
                        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                else
                    [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,tol,...
                        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
                end
                y0b = c0b/(1/gbar/exp(gnow) - phi/2*(pi0b-pibar)^2);
                rn0b = rnpast^rhor*( rnss*(pi0b/pibar)^psi1*(y0b/ystar)^psi2 )^(1-rhor)*exp(rnow);            

                % expectation term (used when peaflag==1)
                fc0n = beta*c0n^(-tau)/(gamma*exp(znow))/pi0n;
                fp0n = beta*phi*c0n^(-tau)*y0n*(pi0n-pibar)*pi0n;
                fc0b = beta*c0b^(-tau)/(gamma*exp(znow))/pi0b;
                fp0b = beta*phi*c0b^(-tau)*y0b*(pi0b-pibar)*pi0b;
                
            elseif (pfmethod==1)
                
                % current period's c and pi (obtained by current period's fc and fp)
                % in the non-ZLB regime
                fc0n = fcvec0n(i);
                fp0n = fpvec0n(i);
                [c0n pi0n] = pf1(fc0n,fp0n,pibar,invnu,tau,phi);
                y0n = c0n/(1/gbar/exp(gnow) - phi/2*(pi0n-pibar)^2);
                rn0n = rnpast^rhor*( rnss*(pi0n/pibar)^psi1*(y0n/ystar)^psi2 )^(1-rhor)*exp(rnow);

                % in the ZLB regime
                fc0b = fcvec0b(i);
                fp0b = fpvec0b(i);
                [c0b pi0b] = pf1(fc0b,fp0b,pibar,invnu,tau,phi);
                y0b = c0b/(1/gbar/exp(gnow) - phi/2*(pi0b-pibar)^2);
                rn0b = rnpast^rhor*( rnss*(pi0b/pibar)^psi1*(y0b/ystar)^psi2 )^(1-rhor)*exp(rnow);            

                % update the expectation term f with interpolation
                fc0n = 0.0;
                fp0n = 0.0;
                fc0b = 0.0;
                fp0b = 0.0;
                for ighe=1:nghe

                    % next period's fc and fp (obtained by interpolation)                
                    gp = rhog*gnow + xmat(ighe,1);
                    zp = rhoz*znow + xmat(ighe,2);
                    rp = xmat(ighe,3);

                    % in the non-ZLB regime (use rn0n)
                    % first assume the ZLB is not binding, and use coeffcn and
                    % coeffpn
                    fcp = makebas4([rn0n gp zp rp],polyd)*coeffcn;
                    fpp = makebas4([rn0n gp zp rp],polyd)*coeffpn;
                    % next period's c and pi (obtained by next period's fc and fp)                                
                    [c1n pi1n] = pf1(fcp,fpp,pibar,invnu,tau,phi);
                    y1n = c1n/(1/gbar/exp(gp) - phi/2*(pi1n-pibar)^2);
                    rn1n = rn0n^rhor*( rnss*(pi1n/pibar)^psi1*(y1n/ystar)^psi2 )^(1-rhor)*exp(rp);
                
                    % then check if the ZLB is violated by rn1n, and use coeffcb and
                    % coeffpb instead
                    if (rn1n<1.0)

                        fcp = makebas4([rn0n gp zp rp],polyd)*coeffcb;
                        fpp = makebas4([rn0n gp zp rp],polyd)*coeffpb;
                        % next period's c and pi (obtained by next period's fc and fp)                                
                        [c1n pi1n] = pf1(fcp,fpp,pibar,invnu,tau,phi);
                        y1n = c1n/(1/gbar/exp(gp) - phi/2*(pi1n-pibar)^2);
                        rn1n = rn0n^rhor*( rnss*(pi1n/pibar)^psi1*(y1n/ystar)^psi2 )^(1-rhor)*exp(rp);
                        
                    end

                    % in the non-ZLB regime (use rn0b)
                    % first assume the ZLB is not binding, and use coeffcn and
                    % coeffpn
                    fcp = makebas4([rn0b gp zp rp],polyd)*coeffcn;
                    fpp = makebas4([rn0b gp zp rp],polyd)*coeffpn;
                    % next period's c and pi (obtained by next period's fc and fp)                                
                    [c1b pi1b] = pf1(fcp,fpp,pibar,invnu,tau,phi);
                    y1b = c1b/(1/gbar/exp(gp) - phi/2*(pi1b-pibar)^2);
                    rn1b = rn0b^rhor*( rnss*(pi1b/pibar)^psi1*(y1b/ystar)^psi2 )^(1-rhor)*exp(rp);
                
                    % then check if the ZLB is violated by rn1n, and use coeffcb and
                    % coeffpb instead
                    if (rn1b<1.0)

                        fcp = makebas4([rn0b gp zp rp],polyd)*coeffcb;
                        fpp = makebas4([rn0b gp zp rp],polyd)*coeffpb;
                        % next period's c and pi (obtained by next period's fc and fp)                                
                        [c1b pi1b] = pf1(fcp,fpp,pibar,invnu,tau,phi);
                        y1b = c1b/(1/gbar/exp(gp) - phi/2*(pi1b-pibar)^2);
                        rn1b = rn0b^rhor*( rnss*(pi1b/pibar)^psi1*(y1b/ystar)^psi2 )^(1-rhor)*exp(rp);
                        
                    end
                    
                    % current period's fc and fp                
                    fcxn = c1n^(-tau)*rn0n/(gamma*exp(zp))/pi1n; % NOTE: max operator is needed?
                    fpxn = c1n^(-tau)*y1n*(pi1n-pibar)*pi1n/y0n;
                    fcxb = c1b^(-tau)*1.0/(gamma*exp(zp))/pi1b;
                    fpxb = c1b^(-tau)*y1b*(pi1b-pibar)*pi1b/y0b;

                    weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
                    fc0n = fc0n + weight*beta*fcxn;
                    fp0n = fp0n + weight*beta*phi*fpxn;
                    fc0b = fc0b + weight*beta*fcxb;
                    fp0b = fp0b + weight*beta*phi*fpxb;

                end
                
            end
            
            cvec1n(i)  = c0n;
            pivec1n(i) = pi0n;
            yvec1n(i)  = y0n;
            rnvec1n(i) = rn0n;
            fcvec1n(i) = fc0n;
            fpvec1n(i) = fp0n;

            cvec1b(i)  = c0b;
            pivec1b(i) = pi0b;
            yvec1b(i)  = y0b;
            rnvec1b(i) = rn0b;
            fcvec1b(i) = fc0b;
            fpvec1b(i) = fp0b;

        end

        % calculate the norm between the old and new policy functions
        diffcn = max(abs(cvec1n-cvec0n));
        diffpn = max(abs(pivec1n-pivec0n));
        diffrn = max(abs(rnvec1n-rnvec0n));
        diffyn = max(abs(yvec1n-yvec0n));
        diffn  = max([diffcn diffpn diffrn diffyn]);

        diffcb = max(abs(cvec1b-cvec0b));
        diffpb = max(abs(pivec1b-pivec0b));
        diffrb = max(abs(rnvec1b-rnvec0b));
        diffyb = max(abs(yvec1b-yvec0b));
        diffb  = max([diffcb diffpb diffrb diffyb]);

        diff = max([diffn diffb]);

        % update the policy functions
        cvec0n  = damp*cvec0n + (1.0-damp)*cvec1n;
        pivec0n = damp*pivec0n + (1.0-damp)*pivec1n;
        rnvec0n = damp*rnvec0n + (1.0-damp)*rnvec1n;
        yvec0n  = damp*yvec0n + (1.0-damp)*yvec1n;
        fcvec0n = damp*fcvec0n + (1.0-damp)*fcvec1n;
        fpvec0n = damp*fpvec0n + (1.0-damp)*fpvec1n;

        cvec0b  = damp*cvec0b + (1.0-damp)*cvec1b;
        pivec0b = damp*pivec0b + (1.0-damp)*pivec1b;
        rnvec0b = damp*rnvec0b + (1.0-damp)*rnvec1b;
        yvec0b  = damp*yvec0b + (1.0-damp)*yvec1b;
        fcvec0b = damp*fcvec0b + (1.0-damp)*fcvec1b;
        fpvec0b = damp*fpvec0b + (1.0-damp)*fpvec1b;

        % counter for iterations
        iter = iter + 1;

        disp([iter diffn diffb]);

    end

    % generate new EDS grid
    disp(' Simulating the ergodic distribution');
    [sim_rn sim_g sim_z sim_r sim_c sim_pi sim_y] = ...
        calcsimzlb(shocks,coefcn,coefpin,coefrnn,coeffcn,coeffpn,coefcb,coefpib,coeffcb,coeffpb,pfmethod,peaflag,polyd,xmat,wmat,tol,...
        pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss,sigmar,sigmag,sigmaz);

    disp(' Generating the EDS grid');
    sim(:,1) = sim_rn;
    sim(:,2) = sim_g;
    sim(:,3) = sim_z;
    sim(:,4) = sim_r;
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
% end of main loop

%% Euler errors
drop = floor(0.05*simT);
simTT = simT + drop;

rnvec = rnss*ones(simTT,1);
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
    
    % euler errors
    fc0 = 0.0;
    fp0 = 0.0;
    for ighe=1:nghe

        gp = rhog*gnow + xmat(ighe,1);
        zp = rhoz*znow + xmat(ighe,2);
        rp = xmat(ighe,3);

        if (pfmethod==0)

            rnp = makebas4([rn0 gp zp rp],polyd)*coefrnn;
            if (rnp>=1.0)    
                cp = makebas4([rn0 gp zp rp],polyd)*coefcn;
                pip = makebas4([rn0 gp zp rp],polyd)*coefpin;
            else
                cp = makebas4([rn0 gp zp rp],polyd)*coefcb;
                pip = makebas4([rn0 gp zp rp],polyd)*coefpib;
            end

        elseif (pfmethod==1)

            rnp = makebas4([rn0 gp zp rp],polyd)*coefrnn;
            if (rnp>=1.0)    
                fcp = makebas4([rn0 gp zp rp],polyd)*coeffcn;
                fpp = makebas4([rn0 gp zp rp],polyd)*coeffpn;
            else
                fcp = makebas4([rn0 gp zp rp],polyd)*coeffcb;
                fpp = makebas4([rn0 gp zp rp],polyd)*coeffpb;
            end
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
prZLB = sum(rnvec<1.0)/simTT*100;

disp([std(log(dyvec)*100) std(log(pivec)*400) std(log(rnvec)*400) prZLB]);
disp([log10(mean(abs(evec))) log10(max(abs(evec))) t iter]);
%disp([log10(mean(abs(evec(:,1)))) log10(max(abs(evec(:,1)))) t]);

end


% TI with no pea: to be done
% function [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,tol,...
%     pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)
% 
% [x0 rc] = csolve('focnk2PV',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,polyd,xmat,wmat,...
%     pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
% c0 = x0(1);
% pi0 = x0(2);
% 
% end


function [c0 pi0] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,np,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2zlbPVf',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,np,xmat,wmat,...
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