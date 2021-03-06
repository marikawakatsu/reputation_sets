#!/usr/bin/env julia

## ReputationSets.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating reputation evolution and tracking
## when individuals interact in different groups/social contexts.

# Outline:
# 1. Initialize population and set memberships.
# 2. Assign initial reputations for each individual and each set.
# Aggregate individuals' set reputations according to strategy.
# 3. Compute fitness by having each individual act on each other individual
# based on their reputation and the social norm.
# With small probability u, they choose the wrong action.
# 4. Keep track of those actions!
# 5. Choose a random individual to update their strategy via a sigmoid function.
# 6. In the next round, use the output of 4 to generate new reputation scores.
# 7. Goto 3 and repeat.

using Revise, StatsBase, PyPlot, Profile

using ReputationSets

b = 1.0 # benefit to cooperating
c = 0.1 # cost to cooperating

N = 4 # population size
M = 1 # number of sets
K = 1

δ = 0.0 # how much do aggregators favor their own set membership?
ϵ = 0.1 # cutoff for aggregators
w = 1.0/N # selection strength
u_s = 0.1/N # mutation rate between strategies
u_p = 0.0/N # error rate in choosing action
u_a = 0.0/N # error rate in assigning reputation

verbose = true

sets = equal_sets(N, M, K)
game = Game(b, c, δ, ϵ, w, u_s, u_p, u_a, "db")
pop = Population(sets, game, [3,4], verbose)

num_gens = 2
total_interactions = 2.0*sum([length(x) for x in sets.set_pairs])

total_cooperation = Float64[]
fitness_means = Float64[]
strategy_freqs = Array{Float64, 1}[]
strat_fitness_means = Array{Float64, 1}[]

evolve!(pop)

for g in 1:num_gens
	println("$g")
	evolve!(pop, 1)
	push!(total_cooperation, sum(pop.prev_actions)/total_interactions)
	push!(fitness_means, mean(pop.fitnesses))
	gen_freqs = zeros(Float64, 5)
	[gen_freqs[pop.strategies[i]+1] += 1.0/pop.sets.N for i in 1:N]
	push!(strategy_freqs, gen_freqs)
	gen_means = zeros(Float64, 5)
	[gen_means[x+1] += mean(pop.fitnesses[pop.strategies .== x]) for x in 0:2]
	push!(strat_fitness_means, gen_means)
end
#
# strategy_freqs_array = zeros(Float64, num_gens, 3)
# fitness_means_array = zeros(Float64, num_gens, 3)
# for g in 1:num_gens
# 	strategy_freqs_array[g,:] = strategy_freqs[g]
# 	fitness_means_array[g,:] = strat_fitness_means[g]
# end
#
# strat_ids = "compartmentalizer", "forgiving", "draconian"
