% nkpf.m
% Main file for computing the policy functions in small NK model without ZLB
% Written by Takeki Sunakawa
% Last updated: May 28 2018

clear all;
rng('default');

% metaparameters
% pfmethod = 1; % =0: TI, =1: future PEA, =2: current PEA
% cflag = 2; % cross terms = 0: smolyak for np=2, =1: tensor, =2: smolyak
% for np=4 (nv=41)
% np = 4; % order of polynomial, 2 or 4
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

disp(' ');
disp(' order of polynomial, np=2');
%nkpf2(0,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
%nkpf2(1,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
%nkpf2(2,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% 
% disp(' ');
% disp(' order of polynomial, np=2 (Smolyak)');
% nkpf2(0,2,0,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(1,2,0,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
 nkpf2(3,2,1,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);

% disp(' ');
% disp(' order of polynomial, np=4 (Smolyak)');
% nkpf2(0,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(1,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);
% nkpf2(2,4,2,ngh,tol,damp,simT,pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz);

function nkpf2(pfmethod,np,cflag,ngh,tol,damp,simT,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,sigmar,sigmag,sigmaz,rnpct,mr,mg,mz)

tic;

%% setups
peaflag = 1; % =1; expectation term as a whole is approximated (only with pfmethod=0)

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
fcss = beta*css^(-tau)*rnss/gamma/piss;
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

%% main loop
diff = 1e+4;
iter = 0;

while ((diff>tol) && (iter<1000))
    
    % fitting polynomials
    % time iteration
    coefc  = bbtinv*cvec0;
    coefpi = bbtinv*pivec0;
    % future or current PEA
    coeffc = bbtinv*fcvec0;
    coeffp = bbtinv*fpvec0;

    for i=1:nv

        rnpast = rngrid(i);
        gnow = ggrid(i);
        znow = zgrid(i);
        rnow = rgrid(i);
        ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow);
        
        % time iteration
        if (pfmethod==0)
        
            % solve nonlinear equations for c and pi
            x0 = [cvec0(i) pivec0(i)]';
            if (peaflag==1)
                [c0 pi0] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,slopecon,np,cflag,xmat,wmat,tol,...
                    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
            else
                [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,slopecon,np,cflag,xmat,wmat,tol,...
                    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
            end

            y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
            rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);            
            % expectation term (used when peaflag==1)
            fc0 = beta*c0^(-tau)/(gamma*exp(znow))/pi0;
            fp0 = beta*phi*c0^(-tau)*y0*(pi0-pibar)*pi0;
                    
        % future PEA
        elseif (pfmethod==1)

            % current period's c and pi (obtained by current period's fc and fp)
            [c0 pi0] = pf1(fcvec0(i),fpvec0(i),pibar,invnu,tau,phi);            
            y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
            rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);

            % update the expectation term f with interpolation
            xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
            % numerical integral with GH quadrature
            % xz is abscissa and wz is weight for eps
            fc0 = 0.0;
            fp0 = 0.0;
            for ighe=1:nghe

                % next period's fc and fp (obtained by interpolation)
                gp = rhog*gnow + xmat(ighe,1);
                zp = rhoz*znow + xmat(ighe,2);
                rp = xmat(ighe,3);
                xg = slopecon(2,1)*gp + slopecon(2,2);
                xzp = slopecon(3,1)*zp + slopecon(3,2);
                xrp = slopecon(4,1)*rp + slopecon(4,2);
                if (np==2)
                    fcp = poly2s([xrnp xg xzp xrp]',cflag)*coeffc;
                    fpp = poly2s([xrnp xg xzp xrp]',cflag)*coeffp;
                elseif (np==4)
                    fcp = poly4s([xrnp,zp xzp xrp]',cflag)*coeffc;
                    fpp = poly4s([xrnp,zp xzp xrp]',cflag)*coeffp;
                end
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
            
        % current PEA
        elseif (pfmethod==2)
            
            % precomputation of integrals
            xpevec = reshape(xpemat(:,i,:),[np ne]);

            % successive approximation for rn and y
            rn0 = rnvec0(i);
            y0 = yvec0(i);
            [c0 pi0] = pf2(rn0,y0,xpevec,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi);
            y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
            rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);

            % update the current fc and fp by using current c, pi and y
            fc0 = beta*c0^(-tau)/(gamma*exp(znow))/pi0;
            fp0 = beta*phi*c0^(-tau)*y0*(pi0-pibar)*pi0;
                        
        elseif (pfmethod==3)
            
            % successive approximation for rn and y
            rn0 = rnvec0(i);
            y0 = yvec0(i);
            [c0 pi0] = pf2gh(rn0,y0,gnow,znow,xmat,wmat,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz);
            y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
            rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);

            % update the current fc and fp by using current c, pi and y
            fc0 = beta*c0^(-tau)/(gamma*exp(znow))/pi0;
            fp0 = beta*phi*c0^(-tau)*y0*(pi0-pibar)*pi0;
                        
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
    
    % counter for iterations
    iter = iter + 1;
    
    disp([iter diffc diffp diffy]);
    
end

t = toc;
% end of main loop

%% Euler errors
drop = floor(0.05*simT);
simTT = simT + drop;

coefc  = bbtinv*cvec0;
coefpi = bbtinv*pivec0;
coefrn = bbtinv*rnvec0;
coefy  = bbtinv*yvec0;
coeffc = bbtinv*fcvec0;
coeffp = bbtinv*fpvec0;

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
    
    % policy function
    xrn = slopecon(1,1)*rnpast + slopecon(1,2);
    xg = slopecon(2,1)*gnow + slopecon(2,2);
    xz = slopecon(3,1)*znow + slopecon(3,2);
    xr = slopecon(4,1)*rnow + slopecon(4,2);

    if (pfmethod==0)

        if (np==2)
            c0  = poly2s([xrn xg xz xr]',cflag)*coefc;
            pi0 = poly2s([xrn xg xz xr]',cflag)*coefpi;
        elseif (np==4)
            c0  = poly4s([xrn xg xz xr]',cflag)*coefc;
            pi0 = poly4s([xrn xg xz xr]',cflag)*coefpi;
        end
        
    elseif (pfmethod==1)

        if (np==2)
            fc0 = poly2s([xrn xg xz xr]',cflag)*coeffc;
            fp0 = poly2s([xrn xg xz xr]',cflag)*coeffp;
        elseif (np==4)
            fc0 = poly4s([xrn xg xz xr]',cflag)*coeffc;
            fp0 = poly4s([xrn xg xz xr]',cflag)*coeffp;
        end
        
        [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi);
        
    elseif (pfmethod==2)
        
        if (np==2)
            y0  = poly2s([xrn xg xz xr]',cflag)*coefy;
            rn0 = poly2s([xrn xg xz xr]',cflag)*coefrn;
        elseif (np==4)
            y0  = poly4s([xrn xg xz xr]',cflag)*coefy;
            rn0 = poly4s([xrn xg xz xr]',cflag)*coefrn;
        end
        
        % precomputation of integrals
        xgpe = prenormchev(slopecon(2,:),np,gnow,rhog,sigmag);
        xzpe = prenormchev(slopecon(3,:),np,znow,rhoz,sigmaz);
        xrpe = prenormchev(slopecon(4,:),np,rnow,0,sigmar);
        xpevec = [xgpe xzpe xrpe];
        [c0 pi0] = pf2(rn0,y0,xpevec,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi);

    elseif (pfmethod==3)
        
        if (np==2)
            y0  = poly2s([xrn xg xz xr]',cflag)*coefy;
            rn0 = poly2s([xrn xg xz xr]',cflag)*coefrn;
        elseif (np==4)
            y0  = poly4s([xrn xg xz xr]',cflag)*coefy;
            rn0 = poly4s([xrn xg xz xr]',cflag)*coefrn;
        end
        
        [c0 pi0] = pf2gh(rn0,y0,gnow,znow,xmat,wmat,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz);
        
    end

    y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
    rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);
    
    % euler errors
    xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
    fc0 = 0.0;
    fp0 = 0.0;
    for ighe=1:nghe

        % NOTE: numerical integration is used for zp for all the methods
        gp = rhog*gnow + xmat(ighe,1);
        zp = rhoz*znow + xmat(ighe,2);
        rp = xmat(ighe,3);
        xgp = slopecon(2,1)*gp + slopecon(2,2);
        xzp = slopecon(3,1)*zp + slopecon(3,2);
        xrp = slopecon(4,1)*rp + slopecon(4,2);

        if (pfmethod==0)

            if (np==2)
                cp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefc;
                pip = poly2s([xrnp xgp xzp xrp]',cflag)*coefpi;
            elseif (np==4)
                cp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefc;
                pip = poly4s([xrnp xgp xzp xrp]',cflag)*coefpi;
            end
            
        elseif (pfmethod==1)

            if (np==2)
                fcp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffc;
                fpp = poly2s([xrnp xgp xzp xrp]',cflag)*coeffp;
            elseif (np==4)
                fcp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffc;
                fpp = poly4s([xrnp xgp xzp xrp]',cflag)*coeffp;
            end

            [cp pip] = pf1(fcp,fpp,pibar,invnu,tau,phi);
            
        elseif (pfmethod==2)

            if (np==2)
                yp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefy;
                rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrn;
            elseif (np==4)
                yp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefy;
                rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrn;
            end
            
            % precomputation of integrals
            xgpe = prenormchev(slopecon(2,:),np,gp,rhog,sigmag);
            xzpe = prenormchev(slopecon(3,:),np,zp,rhoz,sigmaz);
            xrpe = prenormchev(slopecon(4,:),np,rp,0,sigmar);
            xpevec = [xgpe xzpe xrpe];
            [cp pip] = pf2(rnp,yp,xpevec,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi);            

        elseif (pfmethod==3)

            if (np==2)
                yp  = poly2s([xrnp xgp xzp xrp]',cflag)*coefy;
                rnp = poly2s([xrnp xgp xzp xrp]',cflag)*coefrn;
            elseif (np==4)
                yp  = poly4s([xrnp xgp xzp xrp]',cflag)*coefy;
                rnp = poly4s([xrnp xgp xzp xrp]',cflag)*coefrn;
            end
            
            [cp pip] = pf2gh(rnp,yp,gp,zp,xmat,wmat,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz);            
            
        end
        
        yp = cp/(1/gbar/exp(gp) - phi/2*(pip-pibar)^2);
        fcx = cp^(-tau)*rn0/(gamma*exp(zp))/pip;
        fpx = cp^(-tau)*yp*(pip-pibar)*pip;
        weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
        fc0 = fc0 + weight*fcx;
        fp0 = fp0 + weight*fpx;

    end
    
    evec(time+1,1) = 1 - (c0^tau)*beta*fc0;
    LHS2 = ( (1-invnu)+invnu*c0^tau - phi*(pi0-pibar)*(pi0-.5*invnu*(pi0-pibar)) )*c0^(-tau)*y0;
    evec(time+1,2) = LHS2 + beta*phi*fp0;
    
    rnvec(time+1) = rn0;
    yvec(time+1)  = y0;
    cvec(time+1)  = c0;
    pivec(time+1) = pi0;
    gvec(time+1)  = rhog*gnow + sigmag*randn;
    zvec(time+1)  = rhoz*znow + sigmaz*randn;
    rvec(time+1)  = sigmar*randn;
    
    if(mod(time,100)==0); disp([c0 pi0 time]); end;
    
end

rnvec = rnvec(drop+1:simTT);
dyvec = yvec(drop+1:simTT)./yvec(drop:simTT-1);
cvec  = cvec(drop+1:simTT);
pivec = pivec(drop+1:simTT);

disp([std(log(dyvec)*100) std(log(pivec)*400) std(log(rnvec)*400)]);
disp([log10(mean(abs(evec))) log10(max(abs(evec))) t iter]);
%disp([log10(mean(abs(evec(:,1)))) log10(max(abs(evec(:,1)))) t]);

end

 
function [c0 pi0] = pf0(x0,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,slopecon,np,cflag,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coefc,coefpi,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

end


function [c0 pi0] = pf0f(x0,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,slopecon,np,cflag,xmat,wmat,tol,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

[x0 rc] = csolve('focnk2f',x0,[],tol^2,100,rnpast,gnow,znow,rnow,ystar,coeffc,coeffp,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss);
c0 = x0(1);
pi0 = x0(2);

end


function [c0 pi0] = pf1(fc0,fp0,pibar,invnu,tau,phi)

c0 = fc0^(-1/tau);
% solve quadratic equation for pi
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau;
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end


function [c0 pi0] = pf2(rn0,y0,xpevec,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi)

xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
% precomputation of integrals
% xzpe have integral of basis functions T_j(xp) for j = 1,2,... where xp = d0 + d1*zp at each
% grid point i=1,...,nv

% calculate the current f with interpolation
if (np==2)
    fc0 = poly2sprecomp(xrnp,xpevec,cflag)*coeffc;
    fp0 = poly2sprecomp(xrnp,xpevec,cflag)*coeffp;
elseif (np==4)
    fc0 = poly4sprecomp(xrnp,xpevec,cflag)*coeffc;
    fp0 = poly4sprecomp(xrnp,xpevec,cflag)*coeffp;
end

c0 = (rn0*fc0)^(-1/tau);
% solve quadratic equation for pi
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau/y0; % BUG 180331
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end


function [c0 pi0] = pf2gh(rn0,y0,gnow,znow,xmat,wmat,coeffc,coeffp,slopecon,np,cflag,pibar,invnu,tau,phi,rhog,rhoz)

xrnp = slopecon(1,1)*rn0 + slopecon(1,2);
fc0 = 0.0;
fp0 = 0.0;
for ighe=1:size(xmat,1)

    % next period's fc and fp (obtained by interpolation)
    gp = rhog*gnow + xmat(ighe,1);
    zp = rhoz*znow + xmat(ighe,2);
    rp = xmat(ighe,3);
    xg = slopecon(2,1)*gp + slopecon(2,2);
    xzp = slopecon(3,1)*zp + slopecon(3,2);
    xrp = slopecon(4,1)*rp + slopecon(4,2);
    if (np==2)
        fcx = poly2s([xrnp xg xzp xrp]',cflag)*coeffc;
        fpx = poly2s([xrnp xg xzp xrp]',cflag)*coeffp;
    elseif (np==4)
        fcx = poly4s([xrnp zp xzp xrp]',cflag)*coeffc;
        fpx = poly4s([xrnp zp xzp xrp]',cflag)*coeffp;
    end

    weight = wmat(ighe,1)*wmat(ighe,2)*wmat(ighe,3);
    fc0 = fc0 + weight*fcx;
    fp0 = fp0 + weight*fpx;
    
end

c0 = (rn0*fc0)^(-1/tau);
% solve quadratic equation for pi
a0 = .5*phi*pibar^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau/y0; % BUG 180331
a1 = -phi*pibar*(1-invnu);
a2 = -phi*(1-.5*invnu);
pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2;

end