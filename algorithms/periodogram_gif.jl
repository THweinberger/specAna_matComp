# Script that creates a animated figure of the spectrum of the periodogram
# estimator when trying to resolve two nearby frequencies.
# This highlights the resolution limit of the periodogram estimator.
#
# Parameters:
# N: number of measurements
# f1: fixed frequency of the first sinusoid
# f_ticks: varying frequency of the second sinusoid
# x: vector of the signal amplitudes of the two sinusoids

using Plots
include("periodogram.jl")

N = 10
f1 = 0.50
x =[1,1]
f_ticks=range(0.25,0.75,length=300)
f_range = range(0,1,length=200)

A(f,N) = exp.(2*pi*1im*kron(f,(0:N-1)))

anim = @animate for i=1:length(f_ticks)
    y = A([0.5 f_ticks[i]],N)*x
    spec=periodogram(y,f_range)
    plot(f_range,spec,xlim=(0,1), ylim=(0,1),
        xticks = 0:0.1:1, yticks = 0:0.1:1,
        xlabel="frequency",ylabel="|Y(f)|^2",
        title = "Periodogram",label="")
    plot!([f1 f_ticks[i]],seriestype="vline",
        label=["f1 = $(round(f1,digits=2))",
        "f2 = $(round(f_ticks[i],digits=2))"])
end every 1

gif(anim, string(pwd(),"/periodogram_fps30.gif"), fps = 30)
