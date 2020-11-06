using FFTW

function periodogram(y,f_range)
    # The periodogram estimator computes the estimate of the spectrum of a
    # signal based on multiple samples in time. For line spectral estimation,
    # the estimates of the frequencies of the sinusoids can be found by
    # searching for local maxima in the spectrum.
    #
    # Inputs:
    # y: measurements
    # f_range: sample point of the spectrum \in [0,1)
    #
    # Output:
    # spec: estimated spectrum

    y=y[:,1]
    N = length(y)
    spec(f) = 1/N*abs( sum(y[i]*exp(-1im*2*pi*f*i) for i = 1:N ))^2
    spec = spec.(f_range)/maximum(spec.(f_range))
end
