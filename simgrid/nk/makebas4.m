function X1 = makebas4(sf,d)

x = sf(:,1);
y = sf(:,2);
z = sf(:,3);
w = sf(:,4);
n = size(sf,1);

if (d==1)
    
    X1 = [ones(size(sf,1),1) x y z w];
    
elseif (d==2)
    
    T0 = ones(n,1);
    Tx1 = x;
    Tx2 = x.^2;
    Ty1 = y;
    Ty2 = y.^2;
    Tz1 = z;
    Tz2 = z.^2;
    Tw1 = w;
    Tw2 = w.^2;
    
    X1 = [T0 Tx1 Ty1 Tz1 Tw1 Tx2 Ty2 Tz2 Tw2 Tx1.*Ty1 Tx1.*Tz1 Tx1.*Tw1 Ty1.*Tz1 Ty1.*Tw1 Tz1.*Tw1];

end