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

using Revise, StatsBase

using ReputationSets

b = 1.0 # benefit to cooperating
c = 0.5 # cost to cooperating

N = 5 # population size
M = 2 # number of sets

δ = 0.1 # how much do aggregators favor their own set membership?
ϵ = 0.2 # cutoff for aggregators
w = 1.0 # selection strength
μ = 0.1 # mutation rate between strategies
u = 0.1 # error rate in choosing action

sets = Sets(N, M)
game = Game(b, c, δ, ϵ, w, μ, u)
pop = Population(sets, game, true)

num_gens = 10

for g in 1:num_gens
	evolve!(pop)
end
