using LinearAlgebra

function ESPRIT(Y,K,m=0)
    # Estimation of signal parameters via rotational invariance techniques
    # (ESPRIT) algorithm for line spectral estimation from measurements.
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

    #C = X*X'/L

    if m == 0
    m = ceil(Int,N/2)
    end

    C = zeros(m,m)
    for j = 1:L, i = m:N
        C = C + Y[i-m+1:i,j]*Y[i-m+1:i,j]'/(L*N)
    end
    N = m

    evdC = eigen(C)
    p = sortperm(evdC.values,rev=true)[1:K]

    S = evdC.vectors[:,p]

    phi = S[1:N-1,:]\S[2:N,:]

    f_hat = ((angle.(eigvals(phi)))/(2*pi) .+ 1).%1
    #f_hat = real(log.(eigvals(phi))/(1im*2*pi))
end
