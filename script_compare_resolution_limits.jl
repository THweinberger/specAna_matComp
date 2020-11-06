# The purpose of this script is to visually compare the emirical performance
# of various algorithms for line spetral estimation for separating
# two closeby sinusoids.
#
# Parameters:
# N: number of measurements per snapshot
# K: number of sinusoids
# L: number of snapshots
# m: number of sparse measurements which can be used in some of the algorithms;
#    set m = N for full measurements
# SNR: signal-to-noise-ratio [dB] of the AWGN model
# f: location of the frequencies of the sinusoids

include("MUSIC.jl")
include("rootMUSIC.jl")
include("ESPRIT.jl")
include("TLSesprit.jl")
include("MPM.jl")
include("ANM.jl")
include("EMAC.jl")
include("periodogram.jl")
include("structHMirls.jl")
using Plots, Random

# Parameters
N = 10
K = 2
L = 1
m = 8
SNR = 30

#f = rand(K)
f = [0.3; 0.40]

# Create sparse samples at m uniformly random locations \in [N]
omega = sort(randperm(N)[1:m])
a(f) = exp.(2*pi*1im*f*(0:N-1))
A(f) = hcat(a.(f')...)

# Initialize data matrix
Y = Array{Complex{Float64},2}(undef,N,L)

# Draw a random noise sample from a multivariate Gaussian distribution
noise = (randn(N,L)+1im*randn(N,L))/ (sqrt(2)*10^(SNR/10))

# Draw signal coefficients from the uniform distribution
# on the complex unit sphere
S = exp.(2*pi*rand(K,L)*1im)
#S = [1, 1]

Y = A(f)*S + noise

f_range = range(0, 1, length=1000)
spec_periodogram = periodogram(Y,f_range)
(f_MUSIC, specMUSIC) = MUSIC(Y,K,m)

f_ESPRIT = ESPRIT(Y,K)
f_rootMUSIC = rootMUSIC(Y,K)
f_TLSesprit = TLSesprit(Y,K)
#println(f_TLSesprit)

f_MPM = MPM(Y,K)
#f_ANM = ANM(Y[omega,:],N,K,omega,1/SNR)
#f_EMAC = EMAC(Y[omega,:],N,K,omega,1/SNR)
f_structHMirls = structHMirls(Y[omega,1],N,K,omega,1/SNR)

plot(f,seriestype="vline",label="ground truth",
xlims = (0,1),
ylims = (0,1.1),
xticks = 0:0.1:1,
xlabel = "f",
ylabel = "normalized spectra"
)
#plot!(f_range,spec_periodogram,label="Periodogram")
#plot!(f_range,specMUSIC,label="MUSIC")
plot!(f_rootMUSIC,seriestype="vline",label="root MUSIC")
plot!(f_structHMirls,seriestype="vline",label="structHMirls")
#plot!(f_MUSIC,seriestype="vline",label="MUSIC")
#plot!(f_ESPRIT,seriestype="vline",label="ESPRIT")
plot!(f_TLSesprit,seriestype="vline",label="TLSesprit")
#plot!(f_MPM,seriestype="vline",label="MPM")
#plot!(f_ANM,seriestype="vline",label="ANM")
#plot!(f_EMAC,seriestype="vline",label="EMAC")
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/periodogram_vs_MUSIC")
