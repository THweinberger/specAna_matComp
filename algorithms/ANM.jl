using Convex, SCS

function ANM(Y,N,K,omega,sigma)
# Atomic norm minimization for line spectral estimation from measurements
# of a sparse linear array
#
# Inputs:
# Y: Data matrix
# N: ambient dimension \geq M
# K: number of sinusoids
# omega: indices of the observed entries of the data in the ambient space
# sigma: noise variance of the additive white Gaussian noise with mean zero
#
# M = size(Y,1): number of measurements per snapshot
# L = size(Y,2): number of snapshots
#
# Output:
# f_hat: estimate of the frequencies of the K sinusoids in [rad]

(M,L) = size(Y)
lambda = sqrt(M*(L+log(M) + sqrt(2*L*log(M)))*sigma)
if L == 1
lambda = sqrt(M*log(N)*sigma)
end
X = ComplexVariable(L,L)
Z = ComplexVariable(N,L)
T = ComplexVariable(N,N)
#u = ComplexVariable(2*N-1)

constraints = ([X Z'; Z T] in :SDP)
#constraints += [T[i,j] == u[i-j+N] for i = 1:N, j = 1:N]

for i = 2:N, j=i:N
constraints += [T[i,j] == T[1,j-i+1]]
end

objective = lambda/(2*sqrt(N))*(real(tr(X) + tr(T))) + 1/2* sumsquares(abs(vec(Z[omega,:]-Y)))
problem = minimize(objective, constraints)

solve!(problem, SCS.Optimizer(verbose=true))

f_hat = MPM(Z.value,K)
end
