function f = poly2sprecomp(x,pevec,cflag)

ne = size(pevec,2);

f = 1.0;
fx = [x 2*x^2-1];
f = [f fx];

for ie = 1:ne
    
    fx = pevec(:,ie)';
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
        %index = [1 2];
        index = index2(j2,:);
        if (index(1)==1)
            fx = [x 2*x^2-1];
        else
            fx = pevec(:,index(1)-1)';
        end
        fy = pevec(:,index(2)-1)';
        f = [f fy(1)*fx fy(2)*fx];

    end

    % 3
    for j3 = 1:4
        %index = [1 2 3];
        index = index3(j3,:);
        if (index(1)==1)
            fx = [x 2*x^2-1];
        else
            fx = pevec(:,index(1)-1)';
        end
        fy = pevec(:,index(2)-1)';
        fz = pevec(:,index(3)-1)';
        f = [f fz(1)*[fy(1)*fx fy(2)*fx] fz(2)*[fy(1)*fx fy(2)*fx]];

    end

    % 4
    index = [1 2 3 4];
    fx = [x 2*x^2-1];
    fy = pevec(:,index(2)-1)';
    fz = pevec(:,index(3)-1)';
    fw = pevec(:,index(4)-1)';
    f = [f fw(1)*[fz(1)*[fy(1)*fx fy(2)*fx] fz(2)*[fy(1)*fx fy(2)*fx]] ...
        fw(2)*[fz(1)*[fy(1)*fx fy(2)*fx] fz(2)*[fy(1)*fx fy(2)*fx]]];

end