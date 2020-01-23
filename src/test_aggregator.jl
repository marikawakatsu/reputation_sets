#!/usr/bin/env julia

## test_aggregator.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test file for aggregator_reputation() function in ReputationSets.

using Revise

using ReputationSets

println("case 1:")
println("i and j share no sets in common and δ = 1")
println("g_i should be 1, the default value")
δ, q, r, h_i, h_j, verbose = 1.0, 0.0, [0, 0, 0], [false, true, true], [true, false, false], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)

println("case 2:")
println("i and j share one set in common, δ = 0, q = 0")
println("j's reputations are [0, 0, 0]")
println("numerator should be 0")
println("denominator should be 1")
println("g_i should be 1")
δ, q, r, h_i, h_j, verbose = 0.0, 0.0, [0, 0, 0], [true, true, true], [true, false, false], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)

println("case 3:")
println("i and j share one set in common, δ = 0, q = 0")
println("j's reputations are [1, 1, 1]")
println("numerator should be 1")
println("denominator should be 1")
println("g_i should be 1")
δ, q, r, h_i, h_j, verbose = 0.0, 0.0, [1, 1, 1], [true, true, true], [true, false, false], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)

println("case 4:")
println("i is in sets 1 and 2")
println("j is in sets 1, 2, and 3")
println("δ = 0.5, q = 0")
println("j's reputations are [1, 0, 1]")
println("numerator has three terms")
println("k = 1 term: (1-δ(1-1))r = (1 - 0)r = 1")
println("k = 2 term: (1-δ(1-1))r = (1 - 0)r = 0")
println("k = 3 term: (1-δ(1-0))r = (1 - 0.5)r = 0.5")
println("numerator should be 1.5")
println("denominator has three terms")
println("k = 1 term: 1-δ(1-1) = (1 - 0) = 1")
println("k = 2 term: 1-δ(1-1) = (1 - 0) = 1")
println("k = 3 term: 1-δ(1-0) = (1 - 0.5) = 0.5")
println("denominator should be 2.5")
println("g_i should be 1")
δ, q, r, h_i, h_j, verbose = 0.5, 0.0, [1, 0, 1], [true, true, false], [true, true, true], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)

println("case 5:")
println("same as case 4, but q is now 0.9")
println("g_i should be 0")
δ, q, r, h_i, h_j, verbose = 0.5, 0.9, [1, 0, 1], [true, true, false], [true, true, true], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)

println("case 6:")
println("i is in sets 1, 2, and 3")
println("j is in sets 1, 2, and 4")
println("δ = 0.6, q = 0")
println("j's reputations are [1, 0, 0, 1]")
println("numerator has three terms")
println("k = 1 term: (1-δ(1-1))r = (1 - 0)r = 1")
println("k = 2 term: (1-δ(1-1))r = (1 - 0)r = 0")
println("k = 4 term: (1-δ(1-0))r = (1 - 0.6)r = 0.4")
println("numerator should be 1.4")
println("denominator has three terms")
println("k = 1 term: 1-δ(1-1) = (1 - 0) = 1")
println("k = 2 term: 1-δ(1-1) = (1 - 0) = 1")
println("k = 4 term: 1-δ(1-0) = (1 - 0.6) = 0.4")
println("denominator should be 2.4")
println("g_i should be 1")
δ, q, r, h_i, h_j, verbose = 0.6, 0.0, [1, 0, 0, 1], [true, true, true, false], [true, true, false, true], true
rep = aggregator_reputation(δ, q, r, h_i, h_j, verbose)
