function [Z,Zprob] = tauchen(N,mu,rho,sigma,m)


Z     = zeros(N,1);
Zprob = zeros(N,N);
% constant term
c     = (1-rho)*mu;

% constructing grid point
% e.g., m = 2
zmax  = m*sqrt(sigma^2/(1-rho^2));
zmin  = -zmax;
w = (zmax-zmin)/(N-1);

% evenly spaced grid points
Z = linspace(zmin,zmax,N)';
Z = Z + mu;

% transition matrix
for j = 1:N % today's state
    
    for k = 1:N % tomorrow's state
        
        if k == 1 % at lower bound
        
            Zprob(j,k) = cdf_normal((Z(1)-c-rho*Z(j)+w/2)/sigma);
        
        elseif k == N % at upper bound
        
            Zprob(j,k) = 1 - cdf_normal((Z(N)-c-rho*Z(j)-w/2)/sigma);

        else % usual case
            Zprob(j,k) = cdf_normal((Z(k)-c-rho*Z(j)+w/2)/sigma) - ...
                         cdf_normal((Z(k)-c-rho*Z(j)-w/2)/sigma);
        end
        
    end
    
end


function c = cdf_normal(x)
    c = 0.5*erfc(-x/sqrt(2));

