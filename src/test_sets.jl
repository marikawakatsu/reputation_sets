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

using StatsBase, ReputationSets

b = 1.0 # benefit to cooperating
c = 0.5 # cost to cooperating

N = 4 # population size
M = 2 # number of sets

δ = 0 # how much do aggregators favor their own set membership?
ϵ = 0.2 # cutoff for aggregators
w = 1.0 # selection strength
μ = 0.1 # mutation rate between strategies
u = 0.1 # error rate in choosing action

sets = Sets(N, M)
game = Game(b, c, δ, ϵ, w, μ, u)
pop = Population(sets, game)

for g in 1:num_gens
	evolve!(pop)
	# # compute payoffs
	# payoffs = zeros(Int64, N)
	# actions = Array{Int64, 2}[]
	# for k in 1:M
	# 	for (i, j) in set_pairs[k]
	# 		# individuals cooperate with indvs whose reps are good
	# 		# and defect with indvs whose reps are bad
	# 		if rand() < u
	# 			reputation == 1 ? i_action = 2 : i_action = 1
	# 		else
	# 			reputation == 1 ? i_action = 1 : i_action = 2
	# 		end
	# 		payoffs[i] += game[reputations[i]+1, reputations[j]+1]
	# 		payoffs[j] += game[reputations[j]+1, reputations[i]+1]
	# 	end
	# end
	# # choose individual to reproduce
	# invader, invadee = sample(1:N, 2)
	# if rand() < 1.0/(1.0+exp(-w*(payoffs[invader])-payoffs[invadee]))
	# 	strategies[invadee] = strategies[invader]
	# 	if rand() < μ
	# 		strategies[invadee] = rand(filter!(x->x!=strategies[invader], collect(0:2)))
	# 	end
	#
	# # update reputations
	# indvs_to_compare = sample(1:N, 2)
end
