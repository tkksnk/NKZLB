function f = poly2(x,y,cflag)

fx = [x 2*x^2-1];
fy = [y 2*y^2-1];

if (cflag==1)
    f = [1.0 fx fy fx(1)*fy fx(2)*fy];
else
    f = [1.0 fx fy];
end