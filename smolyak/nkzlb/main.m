% nkpfzlb.m
% Main file for computing the policy functions in small NK model with ZLB
% Written by Takeki Sunakawa
% March 27 2018 @ Pier's Cafe Todoroki, Setagaya
% March 30 2018 @ Co-working space Pao Gotanda
% April 3 2018  @ Doutor Coffee Oyamadai, Setagaya
% May 29 2018 @ Mannheim

clear all;

% metaparameters
% pfmethod = 1; % =0: TI, =1: future PEA, =2: current PEA
% cflag = 2; % cross terms = 0: smolyak for np=2, =1: tensor, =2: smolyak
% for np=4 (nv=41)
np = 2; % order of polynomial, 2 or 4
ngh = 3; % number of gh nodes
tol = 1e-6; % tolerance for main loop
damp = 0.7; % dampling parameter (for old policy functions)
simT = 10000; % length of simulation for Euler error

% parameter values
% these parameters are chosen by me
invnu = 6; % 1/nu = elasticity
gyss = 0.2;
gbar = 1.0/(1.0-gyss);
% from Herbst and Schorfheide (2015) p. 98
rrbar = 0.42;
pibar = 1.0+3.30/400;
gamma = 1.0+0.52/100;
beta  = 1.0/(1.0+rrbar/400);

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
rnpct = 0.1;
mr = 2.0;
mg = 2.0;
mz = 2.0;

% disp(' ');
% disp(' order of polynomial, np=2');
% nkpf2(0,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(1,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(2,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);

% disp(' ');
% disp(' order of polynomial, np=2 (Smolyak)');
% nkpf2(0,2,0,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(1,2,0,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
nkpf2(2,2,0,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% 
% disp(' ');
% disp(' order of polynomial, np=4 (Smolyak)');
% nkpf2(0,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(1,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(2,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);


function nkpf2(pfmethod,np,cflag,ngh,tol,damp,simT,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz)

tic;

%% setups
peaflag = 1; % =1; expectation terms are approximated (only with pfmethod=0)

% number of variables
ns = 4;
ne = ns-1;
nghe = ngh^ne;

% number of grid points
if (cflag==1)
    nv = (np+1)^ns; 
elseif (cflag==2) % only for np=4
    nv = 41;
else
    nv = 1 + np*ns;
end

% steady state values
css = (1-1/invnu)^(1/tau);
piss = pibar;
rnss = gamma*pibar/beta;
yss = gbar*css;
fcss = beta*css^(-tau)/gamma/piss;
fpss = beta*phi*css^(-tau)*yss*(piss-pibar)*piss;

% set up grid points
rnmin = (1-rnpct)*rnss;
rnmax = (1+rnpct)*rnss;
gmin = -mg*sigmag/sqrt(1-rhog^2);
gmax = mg*sigmag/sqrt(1-rhog^2);
zmin = -mz*sigmaz/sqrt(1-rhoz^2);
zmax = mz*sigmaz/sqrt(1-rhoz^2);
rmin = -mr*sigmar;
rmax = mr*sigmar;

if (np==2)
    xgrid = makegrid2(ns,nv,cflag);
elseif (np==4)
    xgrid = makegrid4(ns,nv,cflag);
end

for i = 1:nv

    if (np==2)
        bbt(i,:) = poly2s(xgrid(:,i),cflag);
    elseif (np==4)
        bbt(i,:) = poly4s(xgrid(:,i),cflag);
    end
    
end

bbtinv = inv(bbt);

rngrid = (rnmax-rnmin)/2*xgrid(1,:) + (rnmax+rnmin)/2;
ggrid = (gmax-gmin)/2*xgrid(2,:) + (gmax+gmin)/2;
zgrid = (zmax-zmin)/2*xgrid(3,:) + (zmax+zmin)/2;
rgrid = (rmax-rmin)/2*xgrid(4,:) + (rmax+rmin)/2;

slopecon = zeros(ns,2);
slopecon(1,1) = 2/(rnmax-rnmin);
slopecon(1,2) = -(rnmax+rnmin)/(rnmax-rnmin);
slopecon(2,1) = 2/(gmax-gmin);
slopecon(2,2) = -(gmax+gmin)/(gmax-gmin);
slopecon(3,1) = 2/(zmax-zmin);
slopecon(3,2) = -(zmax+zmin)/(zmax-zmin);
slopecon(4,1) = 2/(rmax-rmin);
slopecon(4,2) = -(rmax+rmin)/(rmax-rmin);

% Gaussian-Hermite quadrature
% xz is abscissa and wz is weight for eps
[xg,wg] = qnwnorm(ngh,0,sigmag);
[xz,wz] = qnwnorm(ngh,0,sigmaz);
[xr,wr] = qnwnorm(ngh,0,sigmar);

wmat = zeros(nghe,ne);
xmat = zeros(nghe,ne);
for ir = 1:ngh

    for iz = 1:ngh

        for ig = 1:ngh

            index = ngh^2*(ir-1)+ngh*(iz-1)+ig;
            xmat(index,:) = [xg(ig) xz(iz) xr(ir)];
            wmat(index,:) = [wg(ig) wz(iz) wr(ir)];

        end

    end

end

% precomputation of integrals
if (pfmethod==2)
    
    xgpe = prenormchev(slopecon(2,:),np,ggrid,rhog,sigmag);
    xzpe = prenormchev(slopecon(3,:),np,zgrid,rhoz,sigmaz);
    xrpe = prenormchev(slopecon(4,:),np,rgrid,0,sigmar);

    xpemat = zeros(np,nv,ne);
    xpemat(:,:,1) = xgpe;
    xpemat(:,:,2) = xzpe;
    xpemat(:,:,3) = xrpe;

end

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
fcvec0b = beta*css^(-tau)/gamma/piss*ones(nv,1); % assumes rn0 = 1.0
fcvec1b = zeros(nv,1);
fpvec0b = fpss*ones(nv,1);
fpvec1b = zeros(nv,1);

%% main loop
diff = 1e+4;
iter = 0;

while ((diff>tol) && (iter<1000))
    
    % fitting polynomials
    % time iteration
    coefcn  = bbtinv*cvec0n;
    coefpin = bbtinv*pivec0n;
    coefcb  = bbtinv*cvec0b;
    coefpib = bbtinv*pivec0b;
    coefrnn = bbtinv*rnvec0n;
    % future or current PEA
    coeffcn = bbtinv*fcvec0n;
    coeffpn = bbtinv*fpvec0n;
    coeffcb = bbtinv*fcvec0b;
    coeffpb = bbtinv*fpvec0b;

    for i=1:nv

        rnpast = rngrid(i);
        gnow = ggrid(i);
        znow = zgrid(i);
        rnow = rgrid(i);
        ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);
        
        if (pfmethod==0)
        
            % solve nonlinear equations for c and pi
            % in the non-ZLB regime
            x0n = [cvec0n(i) pivec0n(i)]';
            if (peaflag==1)
                [c0n pi0n prn] = pf0f(x0n,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,0,slopecon,np,cflag,xmat,wmat,tol,...
                    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
            else
                [c0n pi0n prn] = pf0(x0n,rnpast,gnow,znow,rnow,ystar,coefcn,coefpin,coefcb,coefpib,coefrnn,0,slopecon,np,cflag,xmat,wmat,tol,...
                    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
            end
            y0n = c0n/(1/gbar/exp(gnow) - phi/2*(pi0n-pibar)^2);
            rn0n = rnpast^rhor*( rnss*(pi0n/pibar)^psi1*(y0n/ystar)^psi2 )^(1-rhor)*exp(rnow);            

            % in the ZLB regime
            x0b = [cvec0b(i) pivec0b(i)]';
            if (peaflag==1)
                [c0b pi0b prb] = pf0f(x0b,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,1,slopecon,np,cflag,xmat,wmat,tol,...
                    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
            else
                [c0b pi0b prb] = pf0(x0b,rnpast,gnow,znow,rnow,ystar,coefcn,coefpin,coefcb,coefpib,coefrnn,1,slopecon,np,cflag,xmat,wmat,tol,...
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
            xrnn = slopecon(1,1)*rn0n + slopecon(1,2);
            xrnb = slopecon(1,1)*rn0b + slopecon(1,2);
            fc0n = 0.0;
            fp0n = 0.0;
            prn = 0.0;
            fc0b = 0.0;
            fp0b = 0.0;
            prb = 0.0;
            for ighe=1:nghe

                % next period's fc and fp (obtained by interpolation)                
                gp = rhog*gnow + xmat(ighe,1);
                zp = rhoz*znow + xmat(ighe,2);
                rp = xmat(ighe,3);
                xgp = slopecon(2,1)*gp + slopecon(2,2);
                xzp = slopecon(3,1)*zp + slopecon(3,2);
                xrp = slopecon(4,1)*rp + slopecon(4,2);
                
                % in the non-ZLB regime (use xrnn)
                % first assume the ZLB is not binding, and use coeffcn and
                % coeffpn
                if (np==2)
                    fc0 = poly2s([xrnn xgp xzp xrp]',cflag)*coeffcn;
                    fp0 = poly2s([xrnn xgp xzp xrp]',cflag)*coeffpn;
                    prn = 0.0;
                elseif (np==4)
                    fc0 = poly4s([xrnn zp xzp xrp]',cflag)*coeffcn;
                    fp0 = poly4s([xrnn zp xzp xrp]',cflag)*coeffpn;
                    prn = 0.0;
                end
                % next period's c and pi (obtained by next period's fc and fp)                                
                [c1n pi1n] = pf1(fc0,fp0,pibar,invnu,tau,phi);
                y1n = c1n/(1/gbar/exp(gp) - phi/2*(pi1n-pibar)^2);
                rn1n = rn0n^rhor*( rnss*(pi1n/pibar)^psi1*(y1n/ystar)^psi2 )^(1-rhor)*exp(rp);
                
                % then check if the ZLB is violated by rn1n, and use coeffcb and
                % coeffpb instead
                if (rn1n<1.0)

                    if (np==2)
                        fc0 = poly2s([xrnn xgp xzp xrp]',cflag)*coeffcb;
                        fp0 = poly2s([xrnn xgp xzp xrp]',cflag)*coeffpb;
                        prn = 1.0;
                    elseif (np==4)
                        fc0 = poly4s([xrnn zp xzp xrp]',cflag)*coeffcb;
                        fp0 = poly4s([xrnn zp xzp xrp]',cflag)*coeffpb;
                        prn = 1.0;
                    end
                    % next period's c and pi (obtained by next period's fc and fp)                                
                    [c1n pi1n] = pf1(fc0,fp0,pibar,invnu,tau,phi);
                    y1n = c1n/(1/gbar/exp(gp) - phi/2*(pi1n-pibar)^2);
                    rn1n = rn0n^rhor*( rnss*(pi1n/pibar)^psi1*(y1n/ystar)^psi2 )^(1-rhor)*exp(rp);
                    
                end

                % in the ZLB regime (use xrnb)
                % first assume the ZLB is not binding, and use coeffcn and
                % coeffpn
                if (np==2)
                    fc0 = poly2s([xrnb xgp xzp xrp]',cflag)*coeffcn;
                    fp0 = poly2s([xrnb xgp xzp xrp]',cflag)*coeffpn;
                    prbx = 0.0;
                elseif (np==4)
                    fc0 = poly4s([xrnb zp xzp xrp]',cflag)*coeffcn;
                    fp0 = poly4s([xrnb zp xzp xrp]',cflag)*coeffpn;
                    prbx = 0.0;
                end
                % next period's c and pi (obtained by next period's fc and fp)                
                [c1b pi1b] = pf1(fc0,fp0,pibar,invnu,tau,phi);
                y1b = c1b/(1/gbar/exp(gp) - phi/2*(pi1b-pibar)^2);
                rn1b = rn0b^rhor*( rnss*(pi1b/pibar)^psi1*(y1b/ystar)^psi2 )^(1-rhor)*exp(rp);
                
                % then check if the ZLB is violated by rn1b, and use coeffcb and
                % coeffpb instead
                if (rn1b<1.0)

                    if (np==2)
                        fc0 = poly2s([xrnb xgp xzp xrp]',cflag)*coeffcb;
                        fp0 = poly2s([xrnb xgp xzp xrp]',cflag)*coeffpb;
                        prbx = 1.0;
                    elseif (np==4)
                        fc0 = poly4s([xrnb zp xzp xrp]',cflag)*coeffcb;
                        fp0 = poly4s([xrnb zp xzp xrp]',cflag)*coeffpb;
                        prbx = 1.0;
                    end
                    % next period's c and pi (obtained by next period's fc and fp)                                
                    [c1b pi1b] = pf1(fc0,fp0,pibar,invnu,tau,phi);
                    y1b = c1b/(1/gbar/exp(gp) - phi/2*(pi1b-pibar)^2);
                    rn1b = rn0b^rhor*( rnss*(pi1b/pibar)^psi1*(y1b/ystar)^psi2 )^(1-rhor)*exp(rp);
                    
                end

                % current period's fc and fp                
%                fcxn = c1n^(-tau)*max(rn0n,1.0)/(gamma*exp(zp))/pi1n; % max operator is needed?
                fcxn = c1n^(-tau)*rn0n/(gamma*exp(zp))/pi1n; % NOTE: max operator is needed?
                fpxn = c1n^(-tau)*y1n*(pi1n-pibar)*pi1n/y0n;
                fcxb = c1b^(-tau)*1.0/(gamma*exp(zp))/pi1b;
                fpxb = c1b^(-tau)*y1b*(pi1b-pibar)*pi1b/y0b;
                
                weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
                fc0n = fc0n + weight*beta*fcxn;
                fp0n = fp0n + weight*beta*phi*fpxn;
                prn = prn + weight*prn;
                fc0b = fc0b + weight*beta*fcxb;
                fp0b = fp0b + weight*beta*phi*fpxb;
                prb = prb + weight*prbx;

            end
                        
        % current PEA
        elseif (pfmethod==2)
            
            % precomputation of integrals
            xpevec = reshape(xpemat(:,i,:),[np ne]);

            % in the non-ZLB regime
            % successive approximation for rn and y
            rn0n = rnvec0n(i);
            y0n = yvec0n(i);                   
            [c0n pi0n prn] = pf2(rn0n,y0n,xpevec,gnow,znow,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,0,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar);
            y0n = c0n/(1/gbar/exp(gnow) - phi/2*(pi0n-pibar)^2);
            rn0n = rnpast^rhor*( rnss*(pi0n/pibar)^psi1*(y0n/ystar)^psi2 )^(1-rhor)*exp(rnow);                            

            % in the ZLB regime
            % successive approximation for rn and y
            rn0b = rnvec0b(i);
            y0b = yvec0b(i);
            [c0b pi0b prb] = pf2(rn0b,y0b,xpevec,gnow,znow,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,1,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar);            
            y0b = c0b/(1/gbar/exp(gnow) - phi/2*(pi0b-pibar)^2);
            rn0b = rnpast^rhor*( rnss*(pi0b/pibar)^psi1*(y0b/ystar)^psi2 )^(1-rhor)*exp(rnow);

            % update the current fc and fp by using current c, pi and y
            fc0n = beta*c0n^(-tau)/(gamma*exp(znow))/pi0n;
            fp0n = beta*phi*c0n^(-tau)*y0n*(pi0n-pibar)*pi0n;
            fc0b = beta*c0b^(-tau)/(gamma*exp(znow))/pi0b;
            fp0b = beta*phi*c0b^(-tau)*y0b*(pi0b-pibar)*pi0b;
                        
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
        
        prvecn(i) = prn;
        prvecb(i) = prb;
        
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

t = toc;
% end of main loop

%% Euler errors
drop = floor(0.05*simT);
simTT = simT + drop;

coefcn  = bbtinv*cvec0n;
coefpin = bbtinv*pivec0n;
coefrnn = bbtinv*rnvec0n;
coefyn  = bbtinv*yvec0n;
coeffcn = bbtinv*fcvec0n;
coeffpn = bbtinv*fpvec0n;
coefcb  = bbtinv*cvec0b;
coefpib = bbtinv*pivec0b;
coefrnb = bbtinv*rnvec0b;
coefyb  = bbtinv*yvec0b;
coeffcb = bbtinv*fcvec0b;
coeffpb = bbtinv*fpvec0b;

rnvec = zeros(simTT,1);
gvec = zeros(simTT,1);
zvec = zeros(simTT,1);
rvec = zeros(simTT,1);
rnvec(1) = rnss;
rng(0);

for time = 1:simTT
    
    rnpast = rnvec(time);
    gnow = gvec(time);
    znow = zvec(time);
    rnow = rvec(time);
    ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);
    
    % policy function
    xrn = slopecon(1,1)*rnpast + slopecon(1,2);
    xg = slopecon(2,1)*gnow + slopecon(2,2);
    xz = slopecon(3,1)*znow + slopecon(3,2);
    xr = slopecon(4,1)*rnow + slopecon(4,2);

    if (pfmethod==0)

        if (np==2)

            rn0 = poly2s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                c0  = poly2s([xrn xg xz xr]',cflag)*coefcn;
                pi0 = poly2s([xrn xg xz xr]',cflag)*coefpin;
            else
                c0  = poly2s([xrn xg xz xr]',cflag)*coefcb;
                pi0 = poly2s([xrn xg xz xr]',cflag)*coefpib;                
            end

        elseif (np==4)

            rn0 = poly4s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                c0  = poly4s([xrn xg xz xr]',cflag)*coefcn;
                pi0 = poly4s([xrn xg xz xr]',cflag)*coefpin;
            else
                c0  = poly4s([xrn xg xz xr]',cflag)*coefcb;
                pi0 = poly4s([xrn xg xz xr]',cflag)*coefpib;                
            end

        end
        
    elseif (pfmethod==1)

        if (np==2)
            
            rn0 = poly2s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                fc0 = poly2s([xrn xg xz xr]',cflag)*coeffcn;
                fp0 = poly2s([xrn xg xz xr]',cflag)*coeffpn;
            else
                fc0 = poly2s([xrn xg xz xr]',cflag)*coeffcb;
                fp0 = poly2s([xrn xg xz xr]',cflag)*coeffpb;                
            end
            
        elseif (np==4)

            rn0 = poly4s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                fc0 = poly4s([xrn xg xz xr]',cflag)*coeffcn;
                fp0 = poly4s([xrn xg xz xr]',cflag)*coeffpn;
            else
                fc0 = poly4s([xrn xg xz xr]',cflag)*coeffcb;
                fp0 = poly4s([xrn xg xz xr]',cflag)*coeffpb;                
            end
        
        end
        
        [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi);
        
    elseif (pfmethod==2)
        
        if (np==2)            
            
            rn0 = poly2s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                y0  = poly2s([xrn xg xz xr]',cflag)*coefyn;
            else
                y0  = poly2s([xrn xg xz xr]',cflag)*coefyb;
            end
            
        elseif (np==4)

            rn0 = poly4s([xrn xg xz xr]',cflag)*coefrnn;
            if (rn0>=1.0)
                y0  = poly4s([xrn xg xz xr]',cflag)*coefyn;
            else
                y0  = poly4s([xrn xg xz xr]',cflag)*coefyb;
            end
        
        end
        
        % precomputation of integrals
        xgpe = prenormchev(slopecon(2,:),np,gnow,rhog,sigmag);
        xzpe = prenormchev(slopecon(3,:),np,znow,rhoz,sigmaz);
        xrpe = prenormchev(slopecon(4,:),np,rnow,0,sigmar);
        xpevec = [xgpe xzpe xrpe];
        [c0 pi0] = pf2(rn0,y0,xpevec,gnow,znow,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,0,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar);

    end

    y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
    rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);
    
    % euler errors
    xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
    fc0 = 0.0;
    fp0 = 0.0;
    for ighe=1:nghe

        gp = rhog*gnow + xmat(ighe,1);
        zp = rhoz*znow + xmat(ighe,2);
        rp = xmat(ighe,3);
        xgp = slopecon(2,1)*gp + slopecon(2,2);
        xzp = slopecon(3,1)*zp + slopecon(3,2);
        xrp = slopecon(4,1)*rp + slopecon(4,2);

        if (pfmethod==0)

            if (np==2)

                rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    cp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefcn;
                    pip = poly2s([xrnp xgp xzp xrp]',cflag)*coefpin;
                else
                    cp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefcb;
                    pip = poly2s([xrnp xgp xzp xrp]',cflag)*coefpib;                
                end

            elseif (np==4)

                rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    cp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefcn;
                    pip = poly4s([xrnp xgp xzp xrp]',cflag)*coefpin;
                else
                    cp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefcb;
                    pip = poly4s([xrnp xgp xzp xrp]',cflag)*coefpib;                
                end

            end

        elseif (pfmethod==1)

            if (np==2)

                rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    fcp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffcn;
                    fpp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffpn;
                else
                    fcp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffcb;
                    fpp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffpb;                
                end
                
            elseif (np==4)

                rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    fcp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffcn;
                    fpp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffpn;
                else
                    fcp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffcb;
                    fpp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffpb;                
                end
                
            end

            [cp pip] = pf1(fcp,fpp,pibar,invnu,tau,phi);
            
        elseif (pfmethod==2)

            if (np==2)
                
                rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    yp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefyn;
                else
                    yp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefyb;                    
                end
                
            elseif (np==4)
                
                rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrnn;
                if (rnp>=1.0)
                    yp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefyn;
                else
                    yp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefyb;                    
                end
                
            end
            
            % precomputation of integrals
            xgpe = prenormchev(slopecon(2,:),np,gp,rhog,sigmag);
            xzpe = prenormchev(slopecon(3,:),np,zp,rhoz,sigmaz);
            xrpe = prenormchev(slopecon(4,:),np,rp,0,sigmar);
            xpevec = [xgpe xzpe xrpe];
            [cp pip] = pf2(rnp,yp,xpevec,gp,zp,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,0,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar);
            
        end
        
        yp = cp/(1/gbar/exp(gp) - phi/2*(pip-pibar)^2);
        fcx = cp^(-tau)*rn0/(gamma*exp(zp))/pip;
        fpx = cp^(-tau)*yp*(pip-pibar)*pip;
        weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
        fc0 = fc0 + weight*fcx;
        fp0 = fp0 + weight*fpx;

    end
    
    evec(time,1) = 1 - (c0^tau)*beta*fc0;
    LHS2 = ( (1-invnu)+invnu*c0^tau - phi*(pi0-pibar)*(pi0-.5*invnu*(pi0-pibar)) )*c0^(-tau)*y0;
    evec(time,2) = LHS2 + beta*phi*fp0;
    
    cvec(time+1)  = c0;
    pivec(time+1) = pi0;
    rnvec(time+1) = rn0;
    yvec(time+1)  = y0;
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

 
function [c0 pi0 pr] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefcn,coefpin,coefcb,coefpib,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2ZLB',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coefcn,coefpin,coefcb,coefpib,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

[f pr] = focnk2ZLB(x0,rnpast,gnow,znow,rnow,ystar,coefcn,coefpin,coefcb,coefpib,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);

end


function [c0 pi0 pr] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2ZLBf',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

[f pr] = focnk2ZLBf(x0,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);

end


function [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi)

c0 = fc0^(-1/tau);
% solve quadratic equation for pi
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end


function [c0 pi0 pr] = pf2(rn0,y0,xpevec,gnow,znow,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar)

xrn = slopecon(1,1)*rn0 + slopecon(1,2);
% precomputation of integrals
% xzpe have integral of basis functions T_j(xp) for j = 1,2,... where xp = d0 + d1*zp at each
% grid point i=1,...,nv

% pickup 1st order terms in polynomial
index = [4 6 8];
z = 1 - (coefrnn(1) + coefrnn(2)*xrn ...
+ coefrnn(index(1))*(slopecon(2,1)*rhog*gnow+slopecon(2,2)) ...
+ coefrnn(index(2))*(slopecon(3,1)*rhoz*znow+slopecon(3,2)) ...
+ coefrnn(index(3))*slopecon(4,2));
% cdf of e<=z where e follows N(0,s2)
% sums of independently distributed random variables 
% https://onlinecourses.science.psu.edu/stat414/node/172
s2 = (coefrnn(index(1))*slopecon(2,1)*sigmag)^2 + ...
    (coefrnn(index(2))*slopecon(3,1)*sigmaz)^2 + ...
    (coefrnn(index(3))*slopecon(4,1)*sigmar)^2;
pr = normcdf(z,0,s2^.5);

% calculate the current f with interpolation
if (np==2)    
    fc0n = poly2sprecomp(xrn,xpevec,cflag)*coeffcn;
    fp0n = poly2sprecomp(xrn,xpevec,cflag)*coeffpn;
    fc0b = poly2sprecomp(xrn,xpevec,cflag)*coeffcb;
    fp0b = poly2sprecomp(xrn,xpevec,cflag)*coeffpb;
elseif (np==4)
    fc0n = poly4sprecomp(xrn,xpevec,cflag)*coeffcn;
    fp0n = poly4sprecomp(xrn,xpevec,cflag)*coeffpn;
    fc0b = poly4sprecomp(xrn,xpevec,cflag)*coeffcb;
    fp0b = poly4sprecomp(xrn,xpevec,cflag)*coeffpb;
end

% use the probability of ZLB as the weight
fc0 = (1-pr)*fc0n + pr*fc0b;
fp0 = (1-pr)*fp0n + pr*fp0b;

if (ZLBflag==1)
    c0 = fc0^(-1/tau);
else
    c0 = (rn0*fc0)^(-1/tau);
end
% solve quadratic equation for pi
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau/y0;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end


function [c0 pi0 pr] = pf2gh(rn0,y0,gnow,znow,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,pibar,invnu,tau,phi,rhog,rhoz,sigmag,sigmaz,sigmar)
% current PEA, use GH quadrature instead

nghe = size(xmat,1);

xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
fc0 = 0.0;
fp0 = 0.0;
pr = 0.0;
for ighe=1:nghe

    gp = rhog*gnow + xmat(ighe,1);
    zp = rhoz*znow + xmat(ighe,2);
    rp = xmat(ighe,3);
    xgp = slopecon(2,1)*gp + slopecon(2,2);
    xzp = slopecon(3,1)*zp + slopecon(3,2);
    xrp = slopecon(4,1)*rp + slopecon(4,2);

    if (np==2)
        
        rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrnn;
        if (rnp>=1.0)
            fcx = poly2s([xrnp xgp xzp xrp]',cflag)*coeffcn;
            fpx = poly2s([xrnp xgp xzp xrp]',cflag)*coeffpn;
            prx = 0.0;
        else
            fcx = poly2s([xrnp xgp xzp xrp]',cflag)*coeffcb;
            fpx = poly2s([xrnp xgp xzp xrp]',cflag)*coeffpb;
            prx = 1.0;
        end
        
    elseif (np==4)

        rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrnn;
        if (rnp>=1.0)
            fcx = poly4s([xrnp xgp xzp xrp]',cflag)*coeffcn;
            fpx = poly4s([xrnp xgp xzp xrp]',cflag)*coeffpn;
            prx = 0.0;
        else
            fcx = poly4s([xrnp xgp xzp xrp]',cflag)*coeffcb;
            fpx = poly4s([xrnp xgp xzp xrp]',cflag)*coeffpb;
            prx = 1.0;
        end
    
    end
    
    weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
    fc0 = fc0 + weight*fcx;
    fp0 = fp0 + weight*fpx;
    pr = pr + weight*prx;

end

if (ZLBflag==1)
    c0 = fc0^(-1/tau);
else
    c0 = (rn0*fc0)^(-1/tau);
end
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau/y0;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end