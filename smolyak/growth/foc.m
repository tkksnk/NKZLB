function f = foc(c0,yterm,znow,coefc,slopecon,np,cflag,xz,wz,tau,beta,alpha,delta,rhoz)

ngh = size(xz,1);

kp = yterm-c0;
xkp = slopecon(1,1)*kp + slopecon(1,2);
f0 = 0;
for igh=1:ngh

    zp = rhoz*znow + xz(igh);
    xzp = slopecon(2,1)*zp + slopecon(2,2);
    if (np==2)
        cp = poly2(xkp,xzp,cflag)*coefc;
    elseif (np==4)
        cp = poly4(xkp,xzp,cflag)*coefc;
    end
    f0 = f0 + wz(igh)*(alpha*exp(zp)*kp^(alpha-1)+1-delta)*(1/cp)^tau;

end

f = (1/c0)^tau - beta*f0;