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

module ReputationSets

export Sets, Game, Population
export evolve!

struct Sets
	# structure for storing information about set membership

	N::Int64 # number of individuals
	M::Int64 # number of sets
	k::Array{Int64, 1} # size of each set
	set_members::Array{Array{Int64, 1}} # list of members of each set
	set_pairs::Array{Tuple{Int64, Int64}, 1} # list of pairs of individuals in each set
	h::Array{Bool, 2} # N by M boolean array denoting set membership

	# constructor for randomized set membership
	function Sets(
		N::Int64,
		M::Int64
		)
		# populate the sets randomly
		h = rand([false, true], N, M)
		# ensure each individual appears in at least one set
		for i in 1:N
			if all([h[i,k] == false for k in 1:M])
				h[i, rand(1:M)] = true
			end
		end
		# also ensure each set has at least two members
		for k in 1:M
			if sum(h[:,k] == 0)
				h[sample(1:N, 2), k] .= true
			elseif sum(h[:,k] == 1)
				h[rand(filter(x -> h[x,k] != 1, 1:N)), k] = true
			end
		end
		set_members = [Int64[] for k in 1:M]
		for indices in findall(x->x==true, h)
			push!(set_members[indices[2]],indices[1])
		end
		# set_pairs is a list of tuples
		# each tuple is a pair of individuals in that set
		set_pairs = [Tuple{Int64, Int64}[] for k in 1:M]
		for k in 1:M
			for i in set_members[k]
				for j in set_members[k]
					if j < i
						push!(set_pairs[k], (i, j))
					end
				end
			end
		end
		return new(N, M, set_members, set_pairs, h)
	end
end

struct Game
	# structure for storing game parameters and such

	b::Float64 # benefit to cooperating
	c::Float64 # cost to cooperating
	δ::Float64 # how much aggregators favor own-set interactions
	ϵ::Float64 # aggregator cutoff
	w::Float64 # selection strength
	μ::Float64 # mutation rate
	u::Float64 # probability of choosing the "wrong" action
	A::Array{Float64, 2} # the actual game matrix

	function Game(
		b::Float64,
		c::Float64,
		δ::Float64,
		ϵ::Float64,
		w::Float64,
		μ::Float64,
		u::Float64
		)
		return new(b, c, δ, ϵ, w, μ, u, [0 -c; b b-c])
	end
end

mutable struct Population
	# the population of individuals and all information about them
	# this should make it easier to pass around neighbors, etc.

	sets::Sets # not mutable: the set distribution
	game::Game # not mutable: game attributes and parameters
	strategies::Array{Int64, 1} # array of reputation assessment strategies
	reputations::Array{Array{Int64, 1}} # each individual's within-set reputation
	attitudes::Array{Array{Int64, 2}} # each individual's opinion of everyone else in each set
	prev_actions::Array{Array{Int64, 2}} # the last action each individual took toward each other, in each set
	fitnesses::Array{Float64, 1} # array of fitnesses
	generation::Int64
	verbose::Bool # turn this on for error tracking

	# constructor if game and network are already specified
	function Population(
		sets::Sets,
		game::Game
		)
		# begin by initializing the population with random strategies
		strategies = rand(collect(0:2), sets.N)
		# the following attributes are arrays-of-arrays
		# to account for the possibility that sets are not of uniform size
		reputations = [ones(Int64, size(sets.set_members[k])) for k in 1:sets.M]
		attitudes = [ones(Int64, size(sets.set_members[k]), size(sets.set_members[k])) for k in 1:sets.M]
		prev_actions = [ones(Int64, size(sets.set_members[k]), size(sets.set_members[k])) for k in 1:sets.M]
		fitnesses = zeros(Float64, sets.N)
		generation = 0
		return new(network, game, strategies, reputations, fitnesses, generation)
	end
end

function evolve!(
	pop::Population,
	generations::Int64 = 1
	)

	update_fitnesses!(pop)
end

function update_fitnesses!(
	pop::Population
	)

function update_reputations!(
	pop::Population
	)
end

function determine_reputations(
	actions::Array{Array{Int64, 2}},
	set_members::Array{Array{Int64, 1}},
	reputations::Array{Array{Int64, 1}},
	strategies::Array{Int64, 1}
	)
	new_reputations =[zeros(Int64, length(reputations[x])) for x in reputations]
	for k in 1:size(actions)[1]
		for i in 1:size(actions[k])
			j = rand(filter!(x->x!=i, collect(1:size(actions[k]))))
			action = actions[k][i,j]
			new_reputations[k][i] = stern_judging_norm(action, reputations[k][j])
		end
	end
	return new_reputations
end

function determine_attitudes(
	reputations::Array{Array{Int64, 1}},
	strategies::Array{Int64, 1}
	)
	attitudes = zeros(Int64, length(strategies), length(strategies))
	for i in 1:N
		if strategies[i] == 0

		end
	end

end

function aggregator_reputation(
	δ::Float64, # added weight of being in the same set
	q::Float64, # aggregator cutoff
	r::Array{Int64, 1}, # reputation of j in i's eyes, by set
	h_i::Array{Int64, 1}, # set membership of i
	h_j::Array{Int64, 1} # set membership of j
	)
	q = 1-ϵ
	denominator = h_j*(1.0 - δ*(1.0 - h_i))
	numerator = h_j*(1.0 - δ*(1.0 - h_i)) * r
	if numerator/denominator < q
		return 0
	else
		return 1
	end
end

function stern_judging_norm(
	action::Int64, # action taken toward the recipient
	reputation::Int64 # reputation of recipient
	)
	norm = [1 0; 0 1]
	return norm[reputation+1, action+1]
end

end

# final end statement to close the module
