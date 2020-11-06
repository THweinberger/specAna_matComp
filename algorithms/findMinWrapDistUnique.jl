using Combinatorics
include("WrapDist.jl")

function findMinWrapDistUnique(f_true,f_hat)
# Assign each of the ground truth frequencies to closest
# (squared error) canditate of estimated frequencies.
# Metric is the wrap-aroud distance on the torus $\mathds{T}=[0,1]$.
# Note: each frequency in f_hat can only be assigned to the closest frequency
# in f_true and can only be used to explain multiple components
# (see findMinWrapDist).
# Note: Complexity of O(n!) where n=length(f_true)!
#
# Inputs:
# f_true: ground truth frequencies
# f_hat: estimated frequencies
#
# Output:
# dists: wrap-around distances of the frequencies in f_true and f_hat

f_true = f_true[:]
f_hat = f_hat[:]
K = length(f_true)

perms = collect(permutations(f_hat))
dists = Array{Float64,1}(undef,factorial(K))
for ii = 1:factorial(K)
    dists[ii] = norm(WrapDist.(f_true,perms[ii]))
end
(mindist, idx) = findmin(dists)

return WrapDist.(f_true,perms[idx])
end
