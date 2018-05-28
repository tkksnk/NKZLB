function f = poly2precomp(x,ey,cflag)
% ey is precomputed integrals

fx = [x 2*x^2-1];

if (cflag==1)
    f = [1.0 fx ey fx(1)*ey fx(2)*ey];
else
    f = [1.0 fx ey];
end