# The purpose of this script is to compute the empirical mean squared error of
# various algorithms for line spetral estimation for a varying ambient dimension
# of the data while keeping the number of measurements fixed.
# The performance of classical algorithms for data from N measurements is then
# compared to the performance of matrix completion based methods with N sparse
# measurements randomly sampled from the N2 \geq N ambient dimensions.
# This allows to show how much sparse subsampling can increase the performance
# when the number of samples is kept fixed.
#
# Parameters:
# N: number of measurements
# N2range: vector of number of ambient dimensions
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

N2range = ceil.(Int,range(8,32,length=4))
N = 8
L = 1
SNR_db = 60
Niters = 1
f_true = [0.3 0.31]
K = length(f_true)
MSE = zeros(length(N2range),3)

A(f,N) = exp.(2*pi*1im*kron(f,(0:N-1)))
SNR = 10 .^(SNR_db./10)

for ii = 1:length(N2range)
    N2 = N2range[ii]
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

plot(N2range, MSE[:,1],yaxis=:log,tickfontsize = 5,
label="Matrix Pencil Method",
#xlims = (0,1),
#xticks = round.(f_range.-f1; digits=2),
xlabel = "M",
ylabel = "MSE"
)
plot!(N2range,MSE[:,2],label="EMAC")
plot!(N2range,MSE[:,3],label="ANM")
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/mse_N2.png")
