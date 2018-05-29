function xgrid = makegrid2(ns,nv,cflag)

np = 2;
xgrid = zeros(ns,nv);
x = [0 -1 1];

i = 2;

for is=1:ns
    
    for ip=1:np
        
        xgrid(is,i) = x(ip+1);
        i = i+1;
        
    end
    
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
        for i2 = 1:np
            for i1 = 1:np
                xgrid(index(1),i) = x(i1+1);
                xgrid(index(2),i) = x(i2+1);
                i = i+1;
            end
        end

    end

    % 3
    for j3 = 1:4

        index = index3(j3,:);
        for i3 = 1:np
            for i2 = 1:np
                for i1 = 1:np
                    xgrid(index(1),i) = x(i1+1);
                    xgrid(index(2),i) = x(i2+1);
                    xgrid(index(3),i) = x(i3+1);
                    i = i+1;
                end
            end
        end

    end

    % 4
    index = [1 2 3 4];
    for i4 = 1:np
        for i3 = 1:np
            for i2 = 1:np
                for i1 = 1:np
                    xgrid(index(1),i) = x(i1+1);
                    xgrid(index(2),i) = x(i2+1);
                    xgrid(index(3),i) = x(i3+1);
                    xgrid(index(4),i) = x(i4+1);
                    i = i+1;
                end
            end
        end
    end

end