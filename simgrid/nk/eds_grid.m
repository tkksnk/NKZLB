function [X_eds eps] = eds_grid(X0,mu,crit,M)

[simT ns] = size(X0);

% the number of draws
N = floor(simT/mu);
simT = N*mu;
X = zeros(N,ns);

% disp(sprintf(' Selecting each %1d point',mu));

for t = 1:N
    
    X(t,:) = X0(t*mu,:);
    
end

% disp(' done');

muX  = mean(X);
sigX = std(X);
Xn = (X-ones(N,1)*muX)./(ones(N,1)*sigX);

% principle component
%coef = pca(Xn);
[U,S,V] = svd(Xn,0);
coef = V;

PC = Xn*coef;
sigPC = std(PC);
PCn = PC./(ones(N,1)*sigPC);

% disp(' Estimating density function');

%tic
density = calcdens(PCn);
% disp(' ');
%toc

density_sort = sort(density);
% the threshold of density
dbar = density_sort(ceil(crit*N));
index = find(density>dbar);
PCn_sort = PCn(index,:);
% disp(' done');
%PCn_sort = PCn;

% disp(' Constructing grid with bisection');

r1 = min(sqrt(sum(PCn_sort.^2,2)));
r2 = max(sqrt(sum(PCn_sort.^2,2)));
% ehigh = r1/2/M^(1/2);
% elow = r2/(M^(1/2)-1);
ehigh = 1d-4;
elow = 1d+4;
diff = 1d+4;
iter = 0;

while (diff>0 && iter<100)
    
    eps = (ehigh+elow)/2;
    PCn_eds = eds(PCn_sort,eps);

    if (size(PCn_eds,1)>M); ehigh = eps; end;
    if (size(PCn_eds,1)<M); elow  = eps; end;
    diff = abs(size(PCn_eds,1)-M);
    
    iter = iter + 1;
        
end

if (iter==100)
    disp(' The max number of iteration is reached');
    M = size(PCn_eds,1);
end

% disp(' done');

% size(PCn_eds)
% size(sigPC)
PC_eds = PCn_eds.*(ones(M,1)*sigPC);
Xn_eds = PC_eds*inv(coef);
X_eds = Xn_eds.*(ones(M,1)*sigX) + ones(M,1)*muX;


function PCn_eds = eds(PCn,eps)

index = 0;

while(isempty(PCn)==0)
    
    index = index + 1;
    PCn_ = PCn(1,:);
    PCn = PCn(2:end,:);
    dist = sum((ones(size(PCn,1),1)*PCn_-PCn).^2,2).^.5;
    PCn = PCn(find(dist>eps),:);
    PCn_eds(index,:) = PCn_;    
    
end


function density = calcdens(PCn)

N = size(PCn,1);
d = size(PCn,2);

density = zeros(N,1);

hbar = N^(-1/(d+4));
const = 1/(N*(2*pi)^(d/2)*hbar^d);

for t = 1:N
        
%     dist2 = (ones(N,1)*PCn(t,:)-PCn).^2*ones(d,1);
%     s2 = ones(1,N)*exp(-dist2/2/hbar^2);
    dist2 = sum((ones(N,1)*PCn(t,:)-PCn).^2,2);
    s2 = sum(exp(-dist2/2/hbar^2));    
    density(t) = const*s2;
%     if (mod(t,floor(N/40))==0); fprintf('.'); end;
    
end