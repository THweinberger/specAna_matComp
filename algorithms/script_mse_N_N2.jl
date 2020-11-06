# Description: see script_mse_N2 but now we vary N and N2 while keeping
# their ratio fixed s.t. N2 = 1/alpha *N where alpha is the subsampling
# factor.
#
#
# Parameters:
# Nrange: vector of number of measurements
# N2: ambient dimensions
# alpha: subsampling factor
# K: number of sinusoids
# L: number of snapshots
# Niters: number of Monte Carlo iterations
# SNR_db: signal-to-noise-ratio [dB] of the AWGN model
# f_true: location of the frequencies of the sinusoids

include("MUSIC.jl")
include("rootMUSIC.jl")
include("ESPRIT.jl")
include("TLSesprit.jl")
include("MPM.jl")
include("ANM.jl")
include("EMAC.jl")
include("periodogram.jl")
include("findMinWrapDist.jl")
include("findMinWrapDistUnique.jl")
include("WrapDist.jl")

using Plots, Random

Nrange = ceil.(Int,range(6,30,length=5))
L = 1
SNR_db = 30
Niters = 5
alpha = 0.5
f_true = [0.3 0.32 0.34]
K = length(f_true)

MSE = zeros(length(Nrange),3)

A(f,N) = exp.(2*pi*1im*kron(f,(0:N-1)))
SNR = 10 .^(SNR_db./10)

for ii = 1:length(Nrange)
    N = Nrange[ii]
    N2::Int = ceil(1/alpha*N)
    for jj = 1:Niters

        noise = (randn(N,L)+1im*randn(N,L))/ (sqrt(2)*SNR)
        omega = sort(randperm(N2)[1:N])
        S = exp.(2*pi*rand(K,L)*1im)

        Y1 = A(f_true,N)*S + noise
        Y2 = A(f_true,N2)[omega,:]*S + noise

        t0 = time_ns()
        f_MUSIC = MPM(Y1,K)
        t1 = (time_ns()-t0)/1e9
        MinWrapDist = findMinWrapDistUnique(f_true,f_MUSIC)
        MSE[ii,1] += norm(MinWrapDist)^2

        f_EMAC = EMAC(Y2,N2,K,omega,1/SNR)
        MinWrapDist = findMinWrapDistUnique(f_true,f_EMAC)
        MSE[ii,2] += norm(MinWrapDist)^2

        f_ANM = ANM(Y2,N2,K,omega,1/SNR)
        MinWrapDist = findMinWrapDistUnique(f_true,f_ANM)
        MSE[ii,3] += norm(MinWrapDist)^2

    end
end

MSE /= Niters

plot(Nrange, MSE[:,1],yaxis=:log,tickfontsize = 6,
label="Matrix Pencil Method",
#xlims = (0,1),
#xticks = round.(f_range.-f1; digits=2),
xlabel = "N",
ylabel = "MSE"
)
plot!(Nrange,MSE[:,2],label="EMAC")
plot!(Nrange,MSE[:,3],label="ANM")
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/mse_N_N2.png")
