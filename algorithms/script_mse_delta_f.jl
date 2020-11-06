# The purpose of this script is to compute the empirical mean squared error of
# various algorithms for line spetral estimation for separating three closeby
# sinusoids with frequencies [f_1, f_1 + delta, f_1 + 2*delta] for varying
# delta.
#
# Parameters:
# N: number of measurements per snapshot
# L: number of snapshots
# Niters: number of Monte Carlo iterations
# m: number of sparse measurements which can be used in some of the algorithms;
#    set m = N for full measurements
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

N = 10
m = N
Niters = 3
L = 1
SNR_db = 30

SNR = 10^(SNR_db/10)

K = 3
f1 = 0.3
f_range = range(0.0,0.1,length = 5)
MSE = zeros(length(f_range),3)

A(f,N) = exp.(2*pi*1im*kron(f,(0:N-1)))

for ii = 1:length(f_range)
    f_true = [f1 f1+f_range[ii] f1 + 2*f_range[ii]]
    for jj = 1:Niters

        noise = (randn(N,L)+1im*randn(N,L))/ (sqrt(2)*10^(SNR_db/10))
        omega = sort(randperm(m)[1:N])
        S = exp.(2*pi*rand(K,L)*1im)

        Y1 = A(f_true,N)*S + noise
        Y2 = A(f_true,m)[omega,:]*S + noise
        #t0 = time_ns()
        #f_MUSIC = rootMUSIC(Y,K,SNR)
        #t1 = (time_ns()-t0)/1e9
        #MinWrapDist = findMinWrapDistUnique(f_true,f_MUSIC)
        #MSE[ii,1] += norm(MinWrapDist)^2

        t0 = time_ns()
        f_ESPRIT = rootMUSIC(Y1,K)
        t2 = (time_ns()-t0)/1e9
        MinWrapDist = findMinWrapDistUnique(f_true,f_ESPRIT)
        MSE[ii,1] += norm(MinWrapDist)^2

        #t0 = time_ns()
        #f_MPM = MPM(Y,K)
        #t3 = (time_ns()-t0)/1e9
        #MinWrapDist = findMinWrapDistUnique(f_true,f_MPM)
        #MSE[ii,3] += norm(MinWrapDist)^2

        f_EMAC = EMAC(Y2,m,K,omega,1/SNR)
        MinWrapDist = findMinWrapDistUnique(f_true,f_EMAC)
        MSE[ii,2] += norm(MinWrapDist)^2

        f_ANM = ANM(Y2,m,K,omega,1/SNR)
        MinWrapDist = findMinWrapDistUnique(f_true,f_ANM)
        MSE[ii,3] += norm(MinWrapDist)^2

    end
end

MSE /= Niters

plot(f_range, MSE[:,1],yaxis=:log,tickfontsize = 6,
label="root MUSIC",
#xlims = (0,1),
#xticks = round.(f_range.-f1; digits=2),
xlabel = "delta f",
ylabel = "MSE"
)
plot!(f_range,MSE[:,2],label="EMAC")
plot!(f_range,MSE[:,3],label="ANM")
#savefig("/home/thomas/atom-amd64/usr/share/atom/example/mse_delta_f.png")
