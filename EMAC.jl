using Convex, SCS

function EMAC(Y,N,K,omega,sigma)
# Nuclear norm minimization for line spectral estimation from measurements
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
m::Int = ceil(N/2)
n = N + 1 - m
lambda = sigma*L*length(omega)/(1e2*sqrt(N))
#lambda = 1e-4

Q1 = ComplexVariable(n*L,n*L)
Q2 = ComplexVariable(m,m)
H = ComplexVariable(m,n*L)
z = ComplexVariable(N,L)

constraints = ([Q1 H'; H Q2] in :SDP)
constraints += [H[i,j+k*n] == z[i+j-1,k+1] for i = 1:m, j = 1:n, k=0:L-1]
constraints += [sumsquares(abs(vec(z[omega,:]-Y))) < lambda]
objective = real( tr(Q1) + tr(Q2) )
problem = minimize(objective, constraints)

solve!(problem, SCS.Optimizer)

clmnsY0 = setdiff(1:n*L,(1:L).*n)
clmnsY1 = clmnsY0.+1

Y0 = H.value[:,clmnsY0]
Y1 = H.value[:,clmnsY1]

F = svd(Y0)
p = sortperm(F.S,rev=true)[1:K]
f = eigvals(Diagonal(ones(K)./F.S[p])*F.U[:,p]'*Y1*F.V[:,p])
f_hat = (angle.(f)/(2*pi) .+ 1).%1

return f_hat
end
