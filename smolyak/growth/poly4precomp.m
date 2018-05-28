function f = poly4precomp(x,ey,cflag)

fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];

if (cflag==1)
    f = [1.0 fx ey fx(1)*ey fx(2)*ey fx(3)*ey fx(4)*ey]; 
else
    f = [1.0 fx ey];
end