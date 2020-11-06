using LinearAlgebra, PolynomialRoots
include("findlocalmaxima.jl")

function rootMUSIC(Y,K,SNR_db=20,m=0)
    # Estimation of signal parameters via rotational invariance techniques
    # (ESPRIT) algorithm for line spectral estimation from measurements.
    # Can handle multiple snapshots by constructing a low rank Hankel matrix
    # from the data.
    #
    # Inputs:
    # Y: Data matrix
    # K: number of sinusoids
    # SNR_db: SNR in each dimension of the data in dB
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

R = zeros(m,m)
for j = 1:L, i = m:N
    R = R + Y[i-m+1:i,j]*Y[i-m+1:i,j]'/(L*N)
end

    evdR = eigen(R)
    p = sortperm(real(evdR.values))[1:m-K]
    Un = evdR.vectors[:,p]

    M = Un*Un'

    # Construct polynomial
    polcoeffs = Array{Complex{Float64},1}(undef,2*m-1)
    for k = 1:2*m-1
        polcoeffs[k] = sum(diag(M,k-m))
    end

    # Find roots of polynomial
    polyroots = roots(polcoeffs)


    # Remove roots with multiplicity > 1
    idcs = []
    for i = 1:length(polyroots)-1, j = i+1:length(polyroots)
        if abs(polyroots[i]-polyroots[j]) < 0.5/SNR_db + 0.001
            push!(idcs,i)
        end
    end

    polyroots2 = polyroots[setdiff(1:length(polyroots),idcs)]

    # Find K roots closest to unit circle
    dist2UnitCircle = abs.(abs.(polyroots2) .- 1)
    p = sortperm(dist2UnitCircle)[1:K]
    closestRoots = polyroots2[p]

    # Calculate frequencies from roots
    f_hat = (angle.(closestRoots)/(2*pi) .+ 1).%1
end
