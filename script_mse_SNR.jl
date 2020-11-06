# The purpose of this script is to compute the empirical mean squared error of
# various algorithms for line spetral estimation for varying signal-to-noise
# ratios via Monte Carlo iterations.
#
# Parameters:
# N: number of measurements per snapshot
# N2: ambient dimension when N \leq N2 sparse measurements are used
# L: number of snapshots
# Niters: number of Monte Carlo iterations
# SNR_db: vector of signal-to-noise-ratios [dB] of the AWGN model
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

N = 15
N2 = 30
L = 1
SNR_db = range(0,30,length=10)
Niters = 1
f_true = [0.3 0.35 0.40]

K = length(f_true)
MSE = zeros(length(SNR_db),3)

A(f,N) = exp.(2*pi*1im*kron(f,(0:N-1)))
SNR = 10 .^(SNR_db./10)

for jj = 1:Niters
        for ii = 1:length(SNR)
        noise = (randn(N,L)+1im*randn(N,L))/ (sqrt(2)*SNR[ii])
        omega = sort(randperm(N2)[1:N])
        S = exp.(2*pi*rand(K,L)*1im)

        Y1 = A(f_true,N)*S + noise
        Y2 = A(f_true,N2)[omega,:]*S + noise

        t0 = time_ns()
        f_MUSIC = MPM(Y1,K)
        t1 = (time_ns()-t0)/1e9
        MinWrapDist = findMinWrapDistUnique(f_true,f_MUSIC)
        MSE[ii,1] += norm(MinWrapDist)^2

        f_EMAC = EMAC(Y2,N2,K,omega,1/SNR[ii])
        MinWrapDist = findMinWrapDistUnique(f_true,f_EMAC)
        MSE[ii,2] += norm(MinWrapDist)^2

        f_ANM = ANM(Y2,N2,K,omega,1/SNR[ii])
        MinWrapDist = findMinWrapDistUnique(f_true,f_ANM)
        MSE[ii,3] += norm(MinWrapDist)^2

    end
end

MSE /= Niters

plot(SNR_db, MSE[:,1],yaxis=:log,tickfontsize = 6,
label="Matrix Pencil Method",
#xlims = (0,1),
#xticks = round.(f_range.-f1; digits=2),
xlabel = "SNR",
ylabel = "MSE"
)
plot!(SNR_db,MSE[:,2],label="EMAC")
plot!(SNR_db,MSE[:,3],label="ANM")
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/mse_N_N2=30.png")
