function f = poly4s(xx,cflag)

ns = size(xx,1);

f = 1.0;

for is = 1:ns
    
    x = xx(is);
    fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];
    f = [f fx];
    
end

% cross terms
if (cflag==1)

    index2 = [1 2;
        1 3;
        1 4;
        2 3;
        2 4
        3 4];

    index3 = [1 2 3;
        1 2 4;
        1 3 4;
        2 3 4];

    % 2
    for j2 = 1:6

        index = index2(j2,:);
        x = xx(index(1));
        y = xx(index(2));
        fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];
        fy = [y 2*y^2-1 4*y^3-3*y 8*y^4-8*y^2-1];
        f = [f fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx];

    end

    % 3
    for j3 = 1:4

        index = index3(j3,:);
        x = xx(index(1));
        y = xx(index(2));
        z = xx(index(3));
        fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];
        fy = [y 2*y^2-1 4*y^3-3*y 8*y^4-8*y^2-1];
        fz = [z 2*z^2-1 4*z^3-3*z 8*z^4-8*z^2-1];
        f = [f fz(1)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(2)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] ...
               fz(3)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(4)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx]];

    end

    % 4
    index = [1 2 3 4];
    x = xx(index(1));
    y = xx(index(2));
    z = xx(index(3));
    w = xx(index(4));
    fx = [x 2*x^2-1 4*x^3-3*x 8*x^4-8*x^2-1];
    fy = [y 2*y^2-1 4*y^3-3*y 8*y^4-8*y^2-1];
    fz = [z 2*z^2-1 4*z^3-3*z 8*z^4-8*z^2-1];
    fw = [w 2*w^2-1 4*w^3-3*w 8*w^4-8*w^2-1];
    f = [f fw(1)*[fz(1)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(2)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] ...
            fz(3)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(4)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx]] ...
           fw(2)*[fz(1)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(2)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] ...
            fz(3)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(4)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx]] ...
           fw(3)*[fz(1)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(2)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] ...
            fz(3)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(4)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx]] ...
           fw(4)*[fz(1)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(2)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] ...
            fz(3)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx] fz(4)*[fy(1)*fx fy(2)*fx fy(3)*fx fy(4)*fx]]];

elseif (cflag==2)
    
    index2 = [1 2;
        1 3;
        1 4;
        2 3;
        2 4
        3 4];

    % 2
    for j2 = 1:6
        %index = [1 2];
        index = index2(j2,:);
        x = xx(index(1));
        y = xx(index(2));
        fx = [x 2*x^2-1];
        fy = [y 2*y^2-1];
        f = [f fy(1)*fx fy(2)*fx];

    end

end