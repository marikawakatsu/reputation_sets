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

	using Random, StatsBase, Combinatorics

	export Sets, Game, Population
	export random_sets, equal_sets
	export evolve!
	export aggregator_reputation
	export update_strategies_db!, update_strategies_pc!
	export update_reputations_and_attitudes!, update_actions_and_fitnesses!

	struct Sets
		# structure for storing information about set membership
		# currently this is immutable, but we can make it mutable later

		N::Int64 # number of individuals
		M::Int64 # number of sets
		set_members::Array{Array{Int64, 1}, 1} # list of members of each set
		set_pairs::Array{Array{Tuple{Int64, Int64}, 1},1} # list of pairs of individuals in each set
		h::Array{Bool, 2} # N by M boolean array denoting set membership
		# call h[:, k] to get all members of set k
		# or h[i, :] to get all of i's set membership

		# constructors for pre-specified set membership
		function Sets(
			h::Array{Bool, 2}
			)
			N, M = size(h)
			# set_members is a list of lists of individuals in each set
			# set_pairs is a list of tuples
			# each tuple is a pair of individuals in that set
			set_members, set_pairs = set_members_and_pairs(h)
			return new(N, M, set_members, set_pairs, h)
		end
		function Sets(
			N::Int64,
			M::Int64,
			h::Array{Bool, 2}
			)
			set_members, set_pairs = set_members_and_pairs(h)
			return new(N, M, set_members, set_pairs, h)
		end
	end

	# constructor for randomized set membership
	function random_sets(
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
		return Sets(h)
	end

	function equal_sets(
		# constructor for each individual belonging to the same
		# number of (equally sized) sets
		N::Int64,
		M::Int64,
		K::Int64 # the number of sets each individual will belong to
		)
		set_ids = collect(combinations(1:M, K))
		rand_indvs = randperm(N)
		h = zeros(Bool, N, M)
		for (i, ii) in enumerate(rand_indvs)
			si = ii%(length(set_ids))+1
			h[i,set_ids[si]] .= true
		end
		return Sets(h)
	end

	struct Game
		# structure for storing game parameters and such

		b::Float64 # benefit to cooperating
		c::Float64 # cost to cooperating
		δ::Float64 # how much aggregators favor own-set interactions
		ϵ::Float64 # aggregator cutoff
		w::Float64 # selection strength
		u_s::Float64 # mutation rate
		u_p::Float64 # probability of choosing the "wrong" action
		u_a::Float64 # probability of incorrectly assigning reputation
		update_rule::String # update rule: PC or death-birth
		A::Array{Float64, 2} # the actual game matrix

		function Game(
			b::Float64,
			c::Float64,
			δ::Float64,
			ϵ::Float64,
			w::Float64,
			u_s::Float64,
			u_p::Float64,
			u_a::Float64,
			update_rule::String
			)
			return new(b, c, δ, ϵ, w, u_s, u_p, u_a,
				update_rule, [0.0 -c; b b-c])
		end

		function Game(
			b::Float64,
			c::Float64,
			δ::Float64,
			ϵ::Float64,
			w::Float64,
			u_s::Float64,
			u_p::Float64,
			u_a::Float64,
			)
			return new(b, c, δ, ϵ, w, u_s, u_p, u_a,
				"death-birth", [0.0 -c; b b-c])
		end
	end

	mutable struct Population
		# the population of individuals and all information about them

		sets::Sets # not mutable: the set distribution
		game::Game # not mutable: game attributes and parameters
		strategies::Array{Int64, 1} # array of reputation assessment strategies
		reputations::Array{Int64, 2} # each individual's public within-set reputation
		attitudes::Array{Int64, 3} # each individual's opinion of everyone else in each set
		prev_actions::Array{Int64, 3} # the last action each individual took toward each other, in each set
		fitnesses::Array{Float64, 1} # array of fitnesses
		generation::Int64 # current generation
		verbose::Bool # turn this on for error tracking

		# constructor if sets and game are already specified
		function Population(
			sets::Sets,
			game::Game
			)
			# begin by initializing the population with random strategies
			strategies = rand(collect(0:2), sets.N)
			reputations = zeros(Int64, sets.M, sets.N) # everyone's reputation, etc. starts at zero
			attitudes = zeros(Int64, sets.M, sets.N, sets.N)
			prev_actions = zeros(Int64, sets.M, sets.N, sets.N)
			fitnesses = zeros(Float64, sets.N)
			generation = 0
			return new(sets, game, strategies, reputations, attitudes, prev_actions, fitnesses, generation, false)
		end

		# same constructor as above, but allows error tracking to be set
		function Population(
			sets::Sets,
			game::Game,
			verbose::Bool
			)
			# begin by initializing the population with random strategies
			strategies = rand(collect(0:2), sets.N)
			reputations = zeros(Int64, sets.M, sets.N) # everyone's reputation, etc. starts at zero
			attitudes = zeros(Int64, sets.M, sets.N, sets.N)
			prev_actions = zeros(Int64, sets.M, sets.N, sets.N)
			fitnesses = zeros(Float64, sets.N)
			generation = 0
			return new(sets, game, strategies, reputations, attitudes, prev_actions, fitnesses, generation, verbose)
		end
	end

	function set_members_and_pairs(
		h::Array{Bool, 2}
		)
		# return a list of list of set members given a boolean array of memberships
		N, M = size(h)
		set_members = [Int64[] for k in 1:M]
		for indices in findall(x->x==true, h)
			push!(set_members[indices[2]],indices[1])
		end
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
		return set_members, set_pairs
	end

	function update_strategies_db!(
		pop::Population
		)
		# randomly choose someone to die
		invadee = sample(1:pop.sets.N)
		# compute the fitnesses of every other individual in the population
		invasion_fitnesses = 1.0 .- pop.game.w .+ pop.game.w*pop.fitnesses[filter(x->x!=invadee, 1:pop.sets.N)]
		# choose a random other individual in the population, weighted by fitness0
		invader = sample(filter(x->x!=invadee, collect(1:pop.sets.N)), Weights(invasion_fitnesses))
		# the chosen individual's strategy replaces the deceased
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("all fitnesses are $(pop.fitnesses)") end
		if rand() < pop.game.u_s
			# this is where we allow the invadee to mutate
			pop.strategies[invadee] = rand(filter(x->x!=pop.strategies[invader], collect(0:2)))
			if pop.verbose println("mutating $invadee to strategy $(pop.strategies[invadee])") end
		else
			pop.strategies[invadee] = pop.strategies[invader]
			if pop.verbose println("adopting strategy $(pop.strategies[invadee])") end
		end
	end

	function update_strategies_pc!(
		pop::Population
		)
		# chooses a random pair of individuals to compare via a sigmoid function
		# the fitter individual has a chance of invading the less fit one
		# (i.e., forcing them to change strategy)
		invader, invadee = sample(1:pop.sets.N, 2)
		# sigmoid update function
		# sanity check: this should be higher if invader fitness > invadee fitness
		update_function = 1.0/(1.0+exp(-pop.game.w*(pop.fitnesses[invader]-pop.fitnesses[invadee])))
		if pop.verbose println("randomly chosen invader $invader and invadee $invadee") end
		if pop.verbose println("fitnesses are $(pop.fitnesses[invader]) and $(pop.fitnesses[invadee])") end
		if pop.verbose println("update function is $update_function") end
		if rand() < update_function
			if rand() < pop.game.u_s
				# this is where we allow the invadee to mutate
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
		# the main evolution function
		# generations allows us to specify how many generations to let the simulation run for
		# we first need to choose actions and update fitnesses
		if pop.verbose println("initiating generation $(pop.generation)") end
		if pop.verbose println("updating actions and fitnesses") end
		update_actions_and_fitnesses!(pop)
		# then make sure everyone's reputation and attitudes are updated
		if pop.verbose println("updating reputations and attitudes") end
		update_reputations_and_attitudes!(pop)
		# then, finally, select a pair of individuals whose fitnesses we will compare
		if pop.verbose println("evolving, generation $(pop.generation)") end
		if pop.game.update_rule ∈ ["pc", "pairwise_comparison", "im", "imitation"]
			update_strategies_pc!(pop)
		elseif pop.game.update_rule ∈ ["db", "death_birth"]
			update_strategies_db!(pop)
		end
	end


	function update_actions_and_fitnesses!(
		pop::Population
		)
		# choose actions between every pair of individuals in every set,
		# then update their fitnesses accordingly

		# initialize all fitnesses and actions at zero
		new_fitnesses = zeros(Float64, pop.sets.N)
		new_actions = zeros(Int64, pop.sets.M, pop.sets.N, pop.sets.N)
		# for each set
		for k in 1:pop.sets.M
			# for each pair within each set
			for (i, j) in pop.sets.set_pairs[k]
				if pop.verbose println("updating actions of $i and $j in set $k") end
				# check i's attitude toward j within set k
				# with probability 1-u_p, they cooperate if j is good and defect if j is bad
				# with probability u_p, the opposite
				rand() > pop.game.u_p ? i_action = pop.attitudes[k,i,j] : i_action = 1 - pop.attitudes[k,i,j]
				# repeat with j
				rand() > pop.game.u_p ? j_action = pop.attitudes[k,j,i] : j_action = 1 - pop.attitudes[k,j,i]
				if pop.verbose println("$i's attitude toward $j is $(pop.attitudes[k,i,j]), so $i does $i_action") end
				if pop.verbose println("$j's attitude toward $i is $(pop.attitudes[k,j,i]), so $j does $j_action") end
				if pop.verbose println("$i earns a payoff of $(pop.game.A[i_action+1, j_action+1])") end
				if pop.verbose println("$j earns a payoff of $(pop.game.A[j_action+1, i_action+1])") end
				# adjust i and j's fitnesses according to their strategies and the game matrix
				new_fitnesses[i] += pop.game.A[i_action+1, j_action+1] # +1 because julia is 1-indexed
				new_fitnesses[j] += pop.game.A[j_action+1, i_action+1]
				# store their last actions toward each other in this set
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
		# update each individual's public reputation within each set
		# by choosing a random action to observe
		# then adjust each individual's private attitude about every other individual

		# initialize reputations and attitudes at zero
		new_reputations = zeros(Int64, pop.sets.M, pop.sets.N)
		new_attitudes = zeros(Int64, pop.sets.M, pop.sets.N, pop.sets.N)
		# for each set
		for k in 1:pop.sets.M
			# for each individual
			for i in 1:length(pop.sets.set_members[k])
				# check a random other individual j and see what i did to j
				j = rand(filter(x -> x != i, pop.sets.set_members[k]))
				action = pop.prev_actions[k,i,j]
				# apply the stern judging norm to determine i's reputation within set k
				rand() > pop.game.u_a ? new_reputations[k,i] = stern_judging_norm(action, pop.reputations[k,j]) : new_reputations[k,i] = 1 - stern_judging_norm(action, pop.reputations[k,j])
				if pop.verbose
					println("in set $k, $i's behavior toward $j is analyzed")
					println("$i did $action and $j's reputation is $(pop.reputations[k,j])")
					println("in set $k, $i's reputation is updated to be $(new_reputations[k,i])")
				end
			end
		end
		pop.reputations = new_reputations
		# again, for each set
		for k in 1:pop.sets.M
			#for each pair of individuals
			for (i, j) in pop.sets.set_pairs[k]
				# update i's attitude toward j in set k and vice versa
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
		# determine the attitude of i toward j in set k
		# based on their set memberships and i's public reputations

		# 0: compartmentalizer
		# individuals get treated according to their within-set reputation

		# 1: forgiving aggregator
		# if an individual has a good reputation in at least one set,
		# their reputation is good overall

		# 2: draconian aggregator
		# if an individual has a bad reputation in at least one set,
		# their reputation is bad overall
		if pop.verbose println("adjusting $i's attitude toward $j in set $k: $i's strategy is $(pop.strategies[i])") end
		if pop.strategies[i] == 0
			# compartmentalizer
			# i's attitude about j in set k is j's reputation in set k
			if pop.verbose println("$i is a compartmentalizer. attitude toward $j in set $k is $(pop.reputations[k,j])") end
			return pop.reputations[k,j]
		elseif pop.strategies[i] == 1
			# forgiving aggregator
			# get j's reputations
			all_reputations = pop.reputations[:,j]
			# based on j's reputations and i, j's set memberships,
			# determine i's overall feeling about j
			# if j is 1 in at least one set, i's attitude about j should be 1
			new_reputation = aggregator_reputation(pop.game.δ, pop.game.ϵ, all_reputations, pop.sets.h[i,:], pop.sets.h[j,:], pop.verbose)
			if pop.verbose println("$i is a forgiving aggregator. $j's reputations look like $all_reputations. attitude is $new_reputation") end
			return new_reputation
		elseif pop.strategies[i] == 2
			# draconian aggregator
			# get j's reputations
			all_reputations = pop.reputations[:,j]
			# based on j's reputations and i, j's set memberships,
			# determine i's overall feeling about j
			# if j is 0 in at least one set, i's attitude about j should be 0
			new_reputation = aggregator_reputation(pop.game.δ, 1-pop.game.ϵ, all_reputations, pop.sets.h[i,:], pop.sets.h[j,:], pop.verbose)
			if pop.verbose println("$i is a draconian aggregator. $j's reputations look like $all_reputations. attitude is $new_reputation") end
			return new_reputation
		end
	end

	function aggregator_reputation(
		δ::Float64, # added weight of being in the same set
		q::Float64, # aggregator cutoff
		r::Array{Int64, 1}, # reputation of j in i's eyes, by set
		h_i::Array{Bool, 1}, # set membership of i
		h_j::Array{Bool, 1}, # set membership of j
		verbose::Bool = false # verbose output
		)
		if sum(h_i .* h_j) == 0 && δ == 1.0
			default_value = 1
			if verbose println("i and j have no sets in common and δ = $δ") end
			if verbose println("returning default value $default_value") end
			return default_value
		else
			if verbose println("i's membership looks like $h_i") end
			if verbose println("j's membership looks like $h_j") end
			if verbose println("we have δ = $δ") end
			denominator = sum(h_j.*(1.0 .- δ*(1.0 .- h_i)))
			numerator = sum(h_j.*(1.0 .- δ*(1.0 .- h_i)) .* r)
			if verbose println("numerator is $numerator") end
			if verbose println("denominator is $denominator") end
			if numerator/denominator < q
				if verbose println("ratio is < $q so H_q = 0") end
				return 0
			else
				if verbose println("ratio is ≧ $q so H_q = 1") end
				return 1
			end
		end
	end

	function stern_judging_norm(
		action::Int64, # action taken toward the recipient
		reputation::Int64 # reputation of recipient
		)
		# the stern judging norm
		# cooperating with good people and defecting with bad people are good
		# anything else is bad
		norm = [1 0; 0 1]
		return norm[reputation+1, action+1] # +1 because julia is 1-indexed
	end

end
# final end statement to close the module
