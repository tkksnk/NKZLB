function xzpe = prenormchev(slopecon,np,zgrid,rhoz,sigmaz)

nv = size(zgrid,2);

d0 = slopecon(2);
d1 = slopecon(1);

% calculate integral of basis functions T_j(xp) for j = 1,2,... where xp = d0 + d1*zp
xzpe = zeros(np,nv);

for i=1:nv

    znow = zgrid(i);
    zpe1 = rhoz*znow;
    zpe2 = (rhoz*znow)^2 + sigmaz^2;
    xzpe(1,i) = d1*zpe1 + d0;
    xzpe(2,i) = 2*(d1^2*zpe2 + 2*d1*d0*zpe1 + d0^2) - 1;

    if (np==4)

        zpe3 = (rhoz*znow)^3 + 3*rhoz*znow*sigmaz^2;    
        zpe4 = (rhoz*znow)^4 + 6*(rhoz*znow)^2*sigmaz^2 + sigmaz^4;    
        xzpe(3,i) = 4*(d1^3*zpe3 + 3*d1^2*d0*zpe2 + 3*d1*d0^2*zpe1 + d0^3) ...
            -3*(d1*zpe1 + d0);
        xzpe(4,i) = 8*(d1^4*zpe4 + 4*d1^3*d0*zpe3 + 6*d1^2*d0^2*zpe2 + 4*d1*d0^3*zpe1 + d0^4) ...
            -8*(d1^2*zpe2 + 2*d1*d0*zpe1 + d0^2) - 1;

    end

end
