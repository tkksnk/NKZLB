function [f pr] = focnk2ZLBf(x0,rnpast,gnow,znow,rnow,ystar,coeffcn,coeffpn,coeffcb,coeffpb,coefrnn,ZLBflag,slopecon,np,cflag,xmat,wmat,...
    pibar,gamma,beta,invnu,gbar,tau,phi,psi1,psi2,rhor,rhog,rhoz,rnss)

nghe = size(xmat,1);

c0  = x0(1,:);
pi0 = x0(2,:);

y0 = c0/(1/gbar/exp(gnow) - phi/2*(pi0-pibar)^2);
rn0 = rnpast^rhor*( rnss*(pi0/pibar)^psi1*(y0/ystar)^psi2 )^(1-rhor)*exp(rnow);

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
    f(1,:) = -c0^(-tau) + fc0;
else
    f(1,:) = -c0^(-tau) + rn0*fc0;
end
f(2,:) = ( (1-invnu)+invnu*c0^tau - phi*(pi0-pibar)*(pi0-.5*invnu*(pi0-pibar)) )*c0^(-tau)*y0 ...
    + fp0;
