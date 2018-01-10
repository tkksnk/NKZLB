% discretionary policy
clear all;

bet = 0.99;
kap = 0.024;
invlam = 1/0.003;
phi = 2.0;
rstar = -log(bet)*100;

pH = 0.01;
pL = 0.75;
sH = rstar;
sL = rstar-5.0;

% (semi-)analytical
A = [1/invlam 0 kap 0;
    kap 0 -1+bet*(1-pH) bet*pH;
    (1-pL) -1+pL (1-pL) pL;
    0 kap bet*(1-pL) -1+bet*pL];
b = [0;0;-sL;0];
x = A\b;

yH  = x(1);
yL  = x(2);
piH = x(3);
piL = x(4);
% check>0
iH = - yH + (1-pH)*(yH+piH) + pH*(yL+piL) + sH

% policy function iteration
Ns = 2;
yvec0 = [yH; yL];
pvec0 = [piH; piL];
ivec0 = [iH; 0];
% yvec0 = zeros(Ns,1);
% pvec0 = zeros(Ns,1);
% ivec0 = zeros(Ns,1);
yvec1 = zeros(Ns,1);
pvec1 = zeros(Ns,1);
ivec1 = zeros(Ns,1);

Gs = [sH; sL];
Ps = [1-pH pH; 1-pL pL]; 

crit = 1e-10;
diff = 1e+4;

while(diff>crit)

    for is = 1:Ns

        s0 = Gs(is);

        ey = Ps(is,:)*yvec0;
        epi = Ps(is,:)*pvec0;

        p0 = bet*epi/(1+kap^2*invlam);
        y0 = -kap*invlam*p0;
        i0 = ey - y0 + epi + s0;
                
        if (i0<0)

            y0 = ey - (0 - epi - s0);
            p0 = kap*y0 + bet*epi;
            i0 = 0;

        end
        
        yvec1(is,1) = y0;
        pvec1(is,1) = p0;
        ivec1(is,1) = i0;

    end
    
    ydiff = max(abs(yvec1-yvec0));
    pdiff = max(abs(pvec1-pvec0));
    idiff = max(abs(ivec1-ivec0));
    diff = max([ydiff pdiff idiff])

    yvec0 = yvec1;
    pvec0 = pvec1;
    ivec0 = ivec1;
    
end