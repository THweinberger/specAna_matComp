using Combinatorics

function findMinWrapDist(f_true,f_hat)
# Assign each of the ground truth frequencies to closest
# (squared error) canditate of estimated frequencies.
# Metric is the wrap-aroud distance on the torus $\mathds{T}=[0,1]$.
#
# Inputs:
# f_true: ground truth frequencies
# f_hat: estimated frequencies
#
# Output:
# dists: wrap-around distances of the frequencies in f_true and f_hat

K = length(f_true)
WrapDist(x,y) = min(abs(x-y),abs(x-y+1))

dists = Array{Float64,1}(undef,K)
for ii = 1:K
    dists[ii] = minimum(WrapDist.(f_true[ii],f_hat))
end

return dists
end
