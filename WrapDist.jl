function WrapDist(x,y)
    # Computes the wrap-around distance on the torus [0,1).
    #
    # Inputs:
    # x: first frequency
    # y: second frequency
    #
    # Output:
    # wrap around distance between x and y
    
    return min(abs(x-y),abs(x-y+1))
end
