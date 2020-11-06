using Plots;
using Gadfly

W_fejer(omega,N) = 1/N * (sin(omega*N/2)/sin(omega/2))^2

N= 30
omega=range(-pi+1e-3,pi-1e-3,length=10000)
W0 = W_fejer.(1e-13,N)
plot(omega,W_fejer.(omega,N)./W0,yaxis=:log,label="",
xlabel="omega",ylabel="W_B(omega)/W(0)")
savefig("/home/thomas/atom-amd64/usr/share/atom/example/FejerKernel")
