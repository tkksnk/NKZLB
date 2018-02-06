function [y dy] = intf1(xx,yy,x)
% 1-dim linear interpolation
nx = size(xx,1);

ix = 1;
% for x
if (x<=xx(1))
    ix = 1;
elseif (x>=xx(nx))
    ix = nx-1;
else
    jx = 2;
    while (jx<=nx)

        if (x<xx(jx))
            ix=jx-1;
            jx=nx;
        end

    jx=jx+1;
    end
end

etax = (x-xx(ix))/(xx(ix+1)-xx(ix));
detax = 1/(xx(ix+1)-xx(ix));
y = (1-etax)*yy(ix) + etax*yy(ix+1);
dy = -detax*yy(ix) + detax*yy(ix+1);