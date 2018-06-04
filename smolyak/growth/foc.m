function f = foc(c0,yterm,znow,coefc,slopecon,np,cflag,xz,wz,tau,beta,alpha,delta,rhoz)
% given c0, returns the value of R(c)

ngh = size(xz,1);

kp = yterm-c0;
xkp = slopecon(1,1)*kp + slopecon(1,2); % transform k to x in [-1,1]
f0 = 0;
% numerical integration with GH quadrature
for igh=1:ngh % index for each quadrature point

    zp = rhoz*znow + xz(igh); % xz is the value of e' at each quadrature point
    xzp = slopecon(2,1)*zp + slopecon(2,2); % transform z to x in [-1,1]
    % interpolation using the coefficients
    if (np==2)
        % second-order polynominal
        cp = poly2(xkp,xzp,cflag)*coefc;
    elseif (np==4)
        % fourth-order polynominal
        cp = poly4(xkp,xzp,cflag)*coefc;
    end
    % the right hand of the euler equation
    f0 = f0 + wz(igh)*(alpha*exp(zp)*kp^(alpha-1)+1-delta)*(1/cp)^tau;

end

% residual function R(c)
f = (1/c0)^tau - beta*f0;