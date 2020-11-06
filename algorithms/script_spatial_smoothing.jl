# This script analyzes the behavior of spatial smoothing by varying the Hankel
# parameter m which determines the trade off between the number of ambient
# measurements and the number of ambient snapshots of the transformed data
# matrix which is constructed in the respective algorithms.
#
# Inputs:
# N: number of measurements
# K: number of sinusoids
# L: number of snapshots
# SNR: SNR in [dB]
# NIters: number of Monte-Carlo iterations
# m: Hankel parameter
# f: vector of the frequencies of the sinusoids

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

N = 30
K = 3
L = 1
SNR = 30
Niters = 100
m = range(20,stop=N)
f_true = [0.3; 0.35;0.4]

MSE = zeros(length(m))

a(f) = exp.(2*pi*1im*f*(0:N-1))
A(f) = hcat(a.(f')...)

for ii = 1:Niters, jj = 1:length(m)
    noise = (randn(N,L)+1im*randn(N,L))/ (sqrt(2)*10^(SNR/10))

    X = exp.(2*pi*rand(K,L)*1im)

    Y = A(f_true)*X + noise

    f_rootMUSIC = TLSesprit(Y,K,m[jj])

    MinWrapDist = findMinWrapDistUnique(f_true,f_rootMUSIC)
    MSE[jj] += norm(MinWrapDist)^2
end

MSE /= Niters*2

plot(m, MSE,yaxis=:log,tickfontsize = 6,
label = "root MUSIC",
xlabel = "window length",
ylabel = "MSE"
)
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/spatial_smoothing.png")
