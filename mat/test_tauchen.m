clear all;

rho = 0.8;
sigma = 0.01;
m = 1.96;
N = 5;
simT = 10000;

[s,P] = tauchen(N,0,rho,sigma,m)

% generate sequence of u by Markov chain
cumP = cumsum(P');
ivec = zeros(simT,1);
ivec(1) = 3;

for t = 1:simT-1
    
    ivec(t+1) = sum(rand-cumP(ivec(t),:)>=0);
    ivec(t+1) = min(ivec(t+1)+1,N);
    
end

svec = s(ivec);

y = svec(2:end);
X = [ones(simT-1,1) svec(1:end-1)];
inv(X'*X)*X'*y
