clear all;

% 2D
xgrid = [-1 0 1]';
ygrid = [-1 0 1]';
for ix = 1:3

    for iy = 1:3
        
        index = 3*(ix-1)+iy;
        x(index) = xgrid(ix);
        y(index) = ygrid(iy);
        
    end
    
end

xsm = [0 -1 1 0 0];
ysm = [0 0 0 -1 1];

figure;
plot(x,y,'b.','MarkerSize',20);
hold on;
plot(xsm,ysm,'ro','MarkerSize',10);
xticks(xgrid);
yticks(ygrid);
% xticks([-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
% yticks([-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
grid on;
axis square;

% 3D
xgrid = [-1 0 1]';
ygrid = [-1 0 1]';
zgrid = [-1 0 1]';
for ix = 1:3

    for iy = 1:3
        
        for iz = 1:3
        
        index = 9*(ix-1)+3*(iy-1)+iz;
        x(index) = xgrid(ix);
        y(index) = ygrid(iy);
        z(index) = zgrid(iz);
        
        end
        
    end
    
end

xsm = [0 -1 1 0 0 0 0];
ysm = [0 0 0 -1 1 0 0];
zsm = [0 0 0 0 0 -1 1];

figure;
plot3(x,y,z,'b.','MarkerSize',20);
hold on;
plot3(xsm,ysm,zsm,'ro','MarkerSize',10);
xticks(xgrid);
yticks(ygrid);
zticks(zgrid);
% xticks([-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
% yticks([-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
grid on;
axis square;