clear all;

x = linspace(-1,1,101)';

% monominal basis
T0 = 1.0*ones(101,1);
T1 = x;
T2 = x.*T1;
T3 = x.*T2;
T4 = x.*T3;
T5 = x.*T4;

figure;
plot(x,T0);
hold on;
plot(x,T1);
plot(x,T2);
plot(x,T3);
plot(x,T4);
plot(x,T5);
legend('T_0','T_1','T_2','T_3','T_4','T_5','Location','SouthEast')
axis square;

% chebyshev basis
T0 = 1.0*ones(101,1);
T1 = x;
T2 = 2*x.*T1 - T0;
T3 = 2*x.*T2 - T1;
T4 = 2*x.*T3 - T2;
T5 = 2*x.*T4 - T3;

figure;
plot(x,T0);
hold on;
plot(x,T1);
plot(x,T2);
plot(x,T3);
plot(x,T4);
plot(x,T5);
x0 = cos(0*pi/5);
x1 = cos(1*pi/5);
x2 = cos(2*pi/5);
x3 = cos(3*pi/5);
x4 = cos(4*pi/5);
x5 = cos(5*pi/5);
disp([x0 x1 x2 x3 x4 x5]);
plot([x0 x0],[-1 1],'k:');
plot([x1 x1],[-1 1],'k:');
plot([x2 x2],[-1 1],'k:');
plot([x3 x3],[-1 1],'k:');
plot([x4 x4],[-1 1],'k:');
plot([x5 x5],[-1 1],'k:');
% x0 = cos(2*(0+1)*pi/2/5);
% x1 = cos(2*(1+1)*pi/2/5);
% x2 = cos(2*(2+1)*pi/2/5);
% x3 = cos(2*(3+1)*pi/2/5);
% x4 = cos(2*(4+1)*pi/2/5);
% x5 = cos(2*(5+1)*pi/2/5);
% disp([x0 x1 x2 x3 x4 x5]);
% plot([x0 x0],[-1 1],'r:');
% plot([x1 x1],[-1 1],'r:');
% plot([x2 x2],[-1 1],'r:');
% plot([x3 x3],[-1 1],'r:');
% plot([x4 x4],[-1 1],'r:');
% plot([x5 x5],[-1 1],'r:');
% legend('T_0','T_1','T_2','T_3','T_4','T_5','Location','SouthEast')
axis square;