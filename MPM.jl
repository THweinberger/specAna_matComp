using LinearAlgebra

function MPM(Y,K, m=0)
    # Matrix pencil method for line spectral estimation from measurements.
    # Can handle multiple snapshots by constructing a low rank Hankel matrix
    # from the data.
    #
    # Inputs:
    # Y: Data matrix
    # K: number of sinusoids
    # m: Hankel parameter of the Hankel matrix of the data
    #
    # Parameters:
    # N = size(Y,1): number of measurements per snapshot
    # L = size(Y,2): number of snapshots
    #
    # Output:
    # f_hat: estimate of the frequencies of the K sinusoids in [rad]

(N,L) = size(Y)

# If not given from the input, use standard pencil parameter m
if m==0
m = floor(Int,N/2)
end

# Construct horizontally stacked Toeplitz matrix
Y_toep = Array{Complex{Float64},2}(undef,N-m,(m+1)*L)

for j = 1:L, i = 1:m+1
    Y_toep[:,i+(j-1)*(m+1)] = Y[i:N-m-1+i,j]
end

clmnsY0 = setdiff(1:(m+1)*L,(1:L).*(m+1))
clmnsY1 = clmnsY0.+1

Y0 = Y_toep[:,clmnsY0]
Y1 = Y_toep[:,clmnsY1]

F = svd(Y0)
p = sortperm(F.S,rev=true)[1:K]

f = eigvals(Diagonal(ones(K)./F.S[p])*F.U[:,p]'*Y1*F.V[:,p])
f = (angle.(f)/(2*pi) .+ 1).%1
end
