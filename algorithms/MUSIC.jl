using LinearAlgebra
include("findlocalmaxima.jl")

function MUSIC(Y,K,m::Int=0,f_steps=1000)
    # Estimation of signal parameters via rotational invariance techniques
    # (ESPRIT) algorithm for line spectral estimation from measurements.
    # Can handle multiple snapshots by constructing a low rank Hankel matrix
    # from the data.
    #
    # Inputs:
    # Y: Data matrix
    # K: number of sinusoids
    # m: Hankel parameter of the Hankel matrix of the data
    # f_steps: number of data points of the MUSIC pseudo spectrum
    #
    # Parameters:
    # N = size(Y,1): number of measurements per snapshot
    # L = size(Y,2): number of snapshots
    #
    # Output:
    # f_hat: estimate of the frequencies of the K sinusoids in [rad]

(N,L) = size(Y)
f_range = range(0, 1, length=f_steps)

if m == 0
m = ceil(Int,N/2)
end

C = zeros(m,m)
for j = 1:L, i = m:N
    C = C + Y[i-m+1:i,j]*Y[i-m+1:i,j]'/(L*N)
end

    evdC = eigen(C)
    p = sortperm(real(evdC.values))[1:m-K]
    Un = evdC.vectors[:,p]

    a(f) = exp.(2*pi*1im*f*(0:m-1))
    spec(f) = 1/norm(Un'*a(f))^2
    Mspec = spec.(f_range)/maximum(spec.(f_range))
    localmaxima = findlocalmaxima(Mspec)
    q = localmaxima[sortperm(Mspec[localmaxima],rev=true)[1:min(K,length(localmaxima))]]
    f_hat = f_range[q]
    return f_hat #, Mspec

end
