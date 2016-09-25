@inline function nextState(E)
  a = 1.0
  σ = 1.0
  a*E + σ*randn()
end

@inline function sample(threshold, discount, dropProb; maxIterations = 10_000)
    E, L, M = 0.0, 0.0, 0.0

    count = 0
    scale = 1.0

    while (-threshold <= E <= threshold || rand() < dropProb) && count <= maxIterations
        count += 1
        L += scale*E*E
        M += scale

        scale *= discount
        E = nextState(E)
    end

    (L, M)
end


@inline function sample_average(threshold, discount, dropProb; iterations::Int=1000)
    ell, emm = 0.0, 0.0

    for i in 1:iterations
        L, M = sample(threshold, discount, dropProb)
        ell += L
        emm += M
    end
    ell /= iterations
    emm /= iterations
    (ell, emm)
end

@fastmath function sa_costly(cost, discount, dropProb ;
    iterations :: Int     = 1_000,
       initial :: Float64 = 1.0,
        decay1 :: Float64 = 0.9,
	    decay2 :: Float64 = 0.999,
	   epsilon :: Float64 = 1e-8, 
        alpha  :: Float64 = 0.01,
             c :: Float64 = 0.1,
	     debug :: Bool    = false, 
    )

    threshold = initial
    trace     = zeros(iterations)

    moment1 = 0.0
    moment2 = 0.0

    weight1 = decay1
    weight2 = decay2

    @inbounds for k in 1:iterations
        threshold_plus  = threshold + c
        threshold_minus = threshold - c

        L_plus , M_plus  = sample_average(threshold_plus,  discount, dropProb)
        L_minus, M_minus = sample_average(threshold_minus, discount, dropProb)

        # We don't do the cost/(1-beta) correction because
        # it gets cancelled when taking the difference
        C_plus  = (L_plus  + cost)/M_plus
        C_minus = (L_minus + cost)/M_minus

        gradient = (C_plus - C_minus)/2c

        moment1 = decay1 * moment1 + (1 - decay1) * gradient
        moment2 = decay2 * moment2 + (1 - decay2) * gradient^2

        corrected1 = moment1/(1 - weight1)
        corrected2 = moment2/(1 - weight2)

        weight1 *= decay1
        weight2 *= decay2

        threshold_delta = corrected1 / ( sqrt(corrected2) + epsilon) 

        threshold -= alpha * threshold_delta

        # If threshold becomes negative, we set it to zero. However, if threshold is
        # zero, M is also zero, and so C is infinity.  
        # So the set the minimum value of threshold to 0.05

        threshold  = max(threshold, 0.05)

        if debug && mod(k,100) == 0
          @printf("#:%8d, threshold=%0.6f\n", k, threshold)
        end
        
        trace[k] = threshold
    end

    return trace
end

@fastmath function sa_constraint(rate, discount, dropProb ;
    iterations :: Int     = 1_000,
       initial :: Float64 = 1.0,
        alpha  :: Float64 = 0.01,
	     debug :: Bool    = false, 
    )

    threshold = initial
    trace     = zeros(iterations)

    target = 1/(1 - discount - rate)

    @inbounds for k in 1:iterations
        _, M = sample_average(threshold, discount, dropProb)

        rl         = alpha/k
        threshold -= rl*(M - target)

        # If threshold becomes negative, we set it to zero. However, if threshold is
        # zero, M is also zero, and so C is infinity.  
        # So the set the minimum value of threshold to 0.05

        threshold  = max(threshold, 0.05)

        if debug && mod(k,100) == 0
          @printf("#:%8d, threshold=%0.6f\n", k, threshold)
        end
        
        trace[k] = threshold
    end

    return trace
end
