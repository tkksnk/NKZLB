function f = poly4(x,y,cflag)

fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];
fy = [y 2*y^2-1 4*y^3-3*y 8*y^4-8*y^2-1];

if (cflag==1)
    f = [1.0 fx fy fx(1)*fy fx(2)*fy fx(3)*fy fx(4)*fy];
else
    f = [1.0 fx fy];
end