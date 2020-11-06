using LinearAlgebra

function TLSesprit(Y,K,m=0)
    # Estimation of signal parameters via rotational invariance techniques
    # (ESPRIT) algorithm for line spectral estimation from measurements.
    # This is an augmented version of ESPRIT which uses a total least-squares
    # (TLS) step for improved performance.
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

    if m == 0
    m = ceil(Int,N/2)
    end

    C = zeros(m,m)
    for j = 1:L
        for i = m : N
            C = C + Y[i-m+1:i,j]*Y[i-m+1:i,j]'/(L*N);
        end
    end

    N = m
    svdC = svd(C)
    V1 = svdC.U[1:N-1,1:K]
    V2 = svdC.U[2:N,1:K]

    #evdC = eigen(C)
    #p = sortperm(real(evdC.values),rev=true)[1:K]
    #V1 = evdC.vectors[1:N-1,p]
    #V2 = evdC.vectors[2:N,p]

    svdV1V2 = svd(hcat(V1,V2))
    E = svdV1V2.V
    E12 = E[1:K,(K+1):2*K]
    E22 = E[(K+1):2*K,(K+1):2*K]

    phi = -E12/E22

    f_hat = ((angle.(eigvals(phi)))/(2*pi) .+ 1).%1
end
