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

	using Random, StatsBase

	export Sets, Game, Population
	export evolve!
	export aggregator_reputation

	struct Sets
		# structure for storing information about set membership
		# currently this is immutable, but we can make it mutable later

		N::Int64 # number of individuals
		M::Int64 # number of sets
		set_members::Array{Array{Int64, 1}, 1} # list of members of each set
		set_pairs::Array{Array{Tuple{Int64, Int64}, 1},1} # list of pairs of individuals in each set
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
				if sum(h[:,k]) == 0
					h[sample(1:N, 2), k] .= true
				elseif sum(h[:,k]) == 1
					h[rand(filter(x -> h[x,k] != true, 1:N)), k] = true
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
						if i < j
							push!(set_pairs[k], (i, j))
						end
					end
				end
			end
			println(set_pairs)
			println(set_members)
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
			return new(b, c, δ, ϵ, w, μ, u, [0.0 -c; b b-c])
		end
	end

	mutable struct Population
		# the population of individuals and all information about them
		# this should make it easier to pass around neighbors, etc.

		sets::Sets # not mutable: the set distribution
		game::Game # not mutable: game attributes and parameters
		strategies::Array{Int64, 1} # array of reputation assessment strategies
		reputations::Array{Int64, 2} # each individual's within-set reputation
		attitudes::Array{Int64, 3} # each individual's opinion of everyone else in each set
		prev_actions::Array{Int64, 3} # the last action each individual took toward each other, in each set
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
			reputations = zeros(Int64, sets.M, sets.N)
			attitudes = zeros(Int64, sets.M, sets.N, sets.N)
			prev_actions = zeros(Int64, sets.M, sets.N, sets.N)
			fitnesses = zeros(Float64, sets.N)
			generation = 0
			return new(sets, game, strategies, reputations, attitudes, prev_actions, fitnesses, generation, false)
		end

		function Population(
			sets::Sets,
			game::Game,
			verbose::Bool
			)
			# begin by initializing the population with random strategies
			strategies = rand(collect(0:2), sets.N)
			# the following attributes are arrays-of-arrays
			# to account for the possibility that sets are not of uniform size
			reputations = zeros(Int64, sets.M, sets.N)
			attitudes = zeros(Int64, sets.M, sets.N, sets.N)
			prev_actions = zeros(Int64, sets.M, sets.N, sets.N)
			fitnesses = zeros(Float64, sets.N)
			generation = 0
			return new(sets, game, strategies, reputations, attitudes, prev_actions, fitnesses, generation, verbose)
		end
	end

	function death_birth!(
		pop::Population
		)
		invader, invadee = sample(1:pop.sets.N, 2)
		update_function = 1.0/(1.0+exp(-pop.game.w*(pop.fitnesses[invader]-pop.fitnesses[invadee])))
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("update function is $update_function") end
		if rand() < update_function
			if rand() < pop.game.μ
				pop.strategies[invadee] = rand(filter(x->x!=pop.strategies[invader], collect(0:2)))
				if pop.verbose println("mutating $invadee to strategy $(pop.strategies[invadee])") end
			else
				pop.strategies[invadee] = pop.strategies[invader]
				if pop.verbose println("adopting strategy $(pop.strategies[invadee])") end
			end
		end
	end

	function evolve!(
		pop::Population,
		generations::Int64 = 1
		)
		println("$(pop.sets.set_pairs), $(pop.sets.set_members)")
		if pop.verbose println("updating actions and fitnesses") end
		update_actions_and_fitnesses!(pop)
		if pop.verbose println("updating reputations and attitudes") end
		update_reputations_and_attitudes!(pop)
		if pop.verbose println("evolving, generation $(pop.generation)") end
		death_birth!(pop)
	end


	function update_actions_and_fitnesses!(
		pop::Population
		)
		new_fitnesses = zeros(Float64, pop.sets.N)
		new_actions = zeros(Int64, pop.sets.M, pop.sets.N, pop.sets.N)
		#println("$new_fitnesses, $new_actions, $(pop.sets.set_pairs)")
		for k in 1:pop.sets.M
			for (i, j) in pop.sets.set_pairs[k]
				if pop.verbose println("updating actions of $i and $j in set $k") end
				rand() > pop.game.u ? i_action = pop.attitudes[k,i,j] : i_action = 1 - pop.attitudes[k,i,j]
				rand() > pop.game.u ? j_action = pop.attitudes[k,j,i] : j_action = 1 - pop.attitudes[k,j,i]
				if pop.verbose println("$i's attitude toward $j is $(pop.attitudes[k,i,j]), so $i does $i_action") end
				if pop.verbose println("$j's attitude toward $i is $(pop.attitudes[k,j,i]), so $j does $j_action") end
				if pop.verbose println("$i earns a payoff of $(pop.game.A[i_action+1, j_action+1])") end
				if pop.verbose println("$j earns a payoff of $(pop.game.A[j_action+1, i_action+1])") end
				new_fitnesses[i] += pop.game.A[i_action+1, j_action+1]
				new_fitnesses[j] += pop.game.A[j_action+1, i_action+1]
				new_actions[k,i,j] = i_action
				new_actions[k,j,i] = j_action
			end
		end
		if pop.verbose println("new fitnesses look like $new_fitnesses") end
		if pop.verbose println("new actions look like $new_actions") end

		pop.fitnesses = new_fitnesses
		pop.prev_actions = new_actions
	end

	function update_reputations_and_attitudes!(
		pop::Population
		)
		new_reputations = zeros(Int64, pop.sets.M, pop.sets.N)
		new_attitudes = zeros(Int64, pop.sets.M, pop.sets.N, pop.sets.N)
		for k in 1:pop.sets.M
			for i in 1:length(pop.sets.set_members[k])
				println("$i, $(pop.sets.set_members[k])")
				j = rand(filter(x -> x != i, pop.sets.set_members[k]))
				action = pop.prev_actions[k,i,j]
				println("$j, $action, $(pop.reputations[k,j])")
				new_reputations[k,i] = stern_judging_norm(action, pop.reputations[k,j])
				if pop.verbose
					println("in set $k, $i's behavior toward $j is analyzed")
					println("$i did $action and $j's reputation is $(pop.reputations[k,j])")
					println("in set $k, $i's reputation is updated to be $(new_reputations[k,i])")
				end
			end
		end
		pop.reputations = new_reputations
		for k in 1:pop.sets.M
			for (i, j) in pop.sets.set_pairs[k]
				new_attitudes[k,i,j] = determine_attitudes(pop, i, j, k)
				new_attitudes[k,j,i] = determine_attitudes(pop, j, i, k)
			end
		end
		pop.attitudes = new_attitudes
	end

	function determine_attitudes(
		pop::Population,
		i::Int64,
		j::Int64,
		k::Int64
		)
		# 0: compartmentalizer.
		# individuals get treated according to their within-set reputation.
		# 1: forgiving aggregator.
		# if an individual has a good reputation in at least one set,
		# their reputation is good overall.
		# 2: draconian aggregator.
		# if an individual has a bad reputation in at least one set,
		# their reputation is bad overall.
		if pop.verbose println("adjusting $i's attitude toward $j in set $k: $i's strategy is $(pop.strategies[i])") end
		if pop.strategies[i] == 0
			if pop.verbose println("$i is a compartmentalizer. attitude toward $j in set $k is $(pop.reputations[k,j])") end
			return pop.reputations[k,j]
		elseif pop.strategies[i] == 1
			all_reputations = pop.reputations[:,j]
			new_reputation = aggregator_reputation(pop.game.δ, pop.game.ϵ, all_reputations, pop.sets.h[i,:], pop.sets.h[j,:])
			if pop.verbose println("$i is a forgiving aggregator. $j's reputations look like $all_reputations. attitude is $new_reputation") end
			return new_reputation
		elseif pop.strategies[i] == 2
			all_reputations = pop.reputations[:,j]
			new_reputation = aggregator_reputation(pop.game.δ, 1-pop.game.ϵ, all_reputations, pop.sets.h[i,:], pop.sets.h[j,:])
			if pop.verbose println("$i is a draconian aggregator. $j's reputations look like $all_reputations. attitude is $new_reputation") end
			return new_reputation
		end
	end

	function aggregator_reputation(
		δ::Float64, # added weight of being in the same set
		q::Float64, # aggregator cutoff
		r::Array{Int64, 1}, # reputation of j in i's eyes, by set
		h_i::Array{Bool, 1}, # set membership of i
		h_j::Array{Bool, 1} # set membership of j
		)
		if sum(h_i .* h_j) == 0
			return 1
		else
			denominator = sum(h_j.*(1.0 .- δ*(1.0 .- h_i)))
			numerator = sum(h_j.*(1.0 .- δ*(1.0 .- h_i)) .* r)
			if numerator/denominator < q
				return 0
			else
				return 1
			end
		end
	end

	# function determine_reputations(
	# 	actions::Array{Array{Int64, 2}},
	# 	set_members::Array{Array{Int64, 1}},
	# 	reputations::Array{Array{Int64, 1}},
	# 	strategies::Array{Int64, 1}
	# 	)
	# 	new_reputations =[zeros(Int64, length(reputations[x])) for x in reputations]
	# 	for k in 1:size(actions)[1]
	# 		for i in 1:size(actions[k])
	# 			j = rand(filter!(x->x!=i, collect(1:size(actions[k]))))
	# 			action = actions[k][i,j]
	# 			new_reputations[k][i] = stern_judging_norm(action, reputations[k][j])
	# 		end
	# 	end
	# 	return new_reputations
	# end

	function stern_judging_norm(
		action::Int64, # action taken toward the recipient
		reputation::Int64 # reputation of recipient
		)
		norm = [1 0; 0 1]
		return norm[reputation+1, action+1]
	end

end
# final end statement to close the module
