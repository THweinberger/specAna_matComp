function findlocalmaxima(signal)
    # Finds the local maxima of a discrete signal by comparing the
    # value at each point to the neighbouring points.
    #
    # Input:
    # signal: denotes the discrete signal from which we seek the local maxima
    #
    # Output:
    # inds: indices of the local maxima, i.e., signal[inds] are the local maxima
    
   inds = Int[]
   if length(signal)>1
       if signal[end]<signal[1]>signal[2]
           push!(inds,1)
       end
       i = 2
       while i < length(signal)
           if signal[i-1]<signal[i]>signal[i+1]
               push!(inds,i)
               i += 1
           end
           i += 1
       end
       if signal[1]<signal[end]>signal[end-1]
           push!(inds,length(signal))
       end
   end
   inds
 end
