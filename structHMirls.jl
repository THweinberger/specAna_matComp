using Convex, SCS, ToeplitzMatrices, LinearAlgebra

function structHMirls(y,N,k,omega,sigma)
    # Line spectral estimator based on matrix completion via an iteratively
    # reweighted least-squares approach.
    # See: Kümmerle, Verdun '19: "Completion of Structured Low-Rank Matrices via
    # Iteratively Reweighted Least Squares".
    #
    # Inputs:
    # Y: Data matrix
    # N: ambient dimension \geq M
    # k: number of sinusoids
    # omega: indices of the observed entries of the data in the ambient space
    # sigma: noise variance of the additive white Gaussian noise with mean zero
    #
    # M = size(Y,1): number of measurements per snapshot
    # L = size(Y,2): number of snapshots
    #
    # Output:
    # f_hat: estimate of the frequencies of the K sinusoids in [rad]

    maxiters = 100
    #lambda = 1e-3
    lambda = sqrt(N*log(length(omega))*sigma)
    #lambda = 1e-1
    alpha = 0
    d1::Int = ceil(N/2)
    d2 = N + 1 - d1

    phi(y) = y[omega]

    function phiStar(z,N,omega)
        x = zeros(Complex{Float64},N)
        x[omega] = z
    end

    function truncSVD(A,k)
        svdA = svd(A)
        return svdA.U[:,1:k]*Diagonal(svdA.S[1:k])*svdA.Vt[1:k,:]
    end

    function kronSum(A,B)
        d1 = size(A)[1]
        d2 = size(B)[1]
        return kron(Matrix{Float64}(I, d2, d2),A) + kron(B,Matrix{Float64}(I, d1, d1))
    end

# Optimization variables
z = ComplexVariable(N)
H = ComplexVariable(d1,d2)
Hvec = ComplexVariable(d1*d2)

# initialize
z1 = phiStar(y,N,omega)
epsilon = svdvals(Hankel(z1[1:d1],z1[d1:end]))[1]
W = Matrix{Float64}(I, d1*d2, d1*d2)
error = Inf
ii = 0
z_new = Inf*ones(N)

# Define OP
constraints = [H[i,j] == z[i+j-1] for i = 1:d1, j = 1:d2]
constraints += [Hvec == vec(H)]
objective = abs(vec(H)'*W*vec(H)) + 1/(2*lambda)*sumsquares(z[omega]-y)
problem = minimize(objective, constraints)

while error > 1e-3 && ii < maxiters
ii += 1
z_old = z_new
if ii == 1
    problem1 = minimize(epsilon^2*sumsquares(z) + 1/(2*lambda)*sumsquares(z[omega]-y), constraints)
    solve!(problem1, SCS.Optimizer)
else
    objective = abs(2*Hvec'*W*Hvec) + 1/(2*lambda)*sumsquares(z[omega]-y)
    problem = minimize(objective, constraints)
    solve!(problem, SCS.Optimizer)
end
z_new = z.value
epsilon = min(epsilon, norm(z_new - z_old)+alpha^(2*k))
error = norm(z_new-z_old)/norm(z_new)
Hz = Hankel(z_new[1:d1],z_new[d1:end])
#Hz = Array{Complex{Float64},2}(undef,d1,d2)
#for i = 1:d1, j = 1:d2
#Hz[i,j] = z_new[i+j-1]
#end
W = inv(kronSum(truncSVD(Hz*Hz',k),truncSVD(Hz'*Hz,k)) + (epsilon^2)*Matrix{Float64}(I, d1*d2, d1*d2))
#println(error)
end

# Retrieve Frequencies
L = 1
clmnsY0 = setdiff(1:d2*L,(1:L).*d2)
clmnsY1 = clmnsY0.+1

Y0 = H.value[:,clmnsY0]
Y1 = H.value[:,clmnsY1]

F = svd(Y0)
p = sortperm(F.S,rev=true)[1:k]

f = eigvals(Diagonal(ones(K)./F.S[p])*F.U[:,p]'*Y1*F.V[:,p])
f_hat = (angle.(f)/(2*pi) .+ 1).%1
#println(f_hat)
#println(svdvals(H.value))
return f_hat

end
