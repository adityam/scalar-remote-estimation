@everywhere begin
    include("packet-drops.jl")

    cost = 100.0
    rate = 0.5

    discount = 0.9
    dropProb = 0.1

    iterations = 10_000

    costly()     =     sa_costly(cost, discount, dropProb; iterations=iterations)
    constraint() = sa_constraint(cost, discount, dropProb; iterations=iterations)
end

function runs(sa; numRuns = 100)
    tuples = pmap(_ -> sa(), 1:numRuns)
    traces = zeros(iterations, numRuns)

    for run in 1:numRuns
        traces[:, run] = tuples[run]
    end

    traces
end

# Use on REPL as follows
# 
# runs(costly; numRuns=100)
# runs(constraint; numRuns=100)
