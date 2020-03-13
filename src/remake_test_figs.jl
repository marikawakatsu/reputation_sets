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

paramcombs = [
(1.0, 0.1, 100, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 3, 2, 0.01, 0.01, 0.1, 0.1),
(1.0, 0.1, 100, 3, 2, 0.01, 0.01, 0.001, 0.001),
(1.0, 0.5, 100, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.9, 100, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 2, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 4, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 100, 4, 3, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 50, 3, 2, 0.01, 0.01, 0.01, 0.01),
(1.0, 0.1, 200, 3, 2, 0.01, 0.01, 0.01, 0.01),
]

for (pci, pc) in enumerate(paramcombs)

	b, c, N, M, K, w, u_s, u_p, u_a = pc

	δ = 0.0
	ϵ = 0.1

	sets = equal_sets(N, M, K)
	game = Game(b, c, δ, ϵ, w, u_s, u_p, u_a, "db")
	pop = Population(sets, game, false)

	num_gens = 10000
	total_interactions = 2.0*sum([length(x) for x in sets.set_pairs])

	total_cooperation = Float64[]
	fitness_means = Float64[]
	strategy_freqs = Array{Float64, 1}[]
	strat_fitness_means = Array{Float64, 1}[]

	evolve!(pop)

	for g in 1:num_gens
		if (g%100 == 0) println("$pci, $g") end
		evolve!(pop, 1)
		push!(total_cooperation, sum(pop.prev_actions)/total_interactions)
		push!(fitness_means, mean(pop.fitnesses))
		gen_freqs = zeros(Float64, 3)
		[gen_freqs[pop.strategies[i]+1] += 1.0/pop.sets.N for i in 1:N]
		push!(strategy_freqs, gen_freqs)
		gen_means = zeros(Float64, 3)
		[gen_means[x+1] += mean(pop.fitnesses[pop.strategies .== x]) for x in 0:2]
		push!(strat_fitness_means, gen_means)
	end

	strategy_freqs_array = zeros(Float64, num_gens, 3)
	fitness_means_array = zeros(Float64, num_gens, 3)
	for g in 1:num_gens
		strategy_freqs_array[g,:] = strategy_freqs[g]
		fitness_means_array[g,:] = strat_fitness_means[g]
	end

	strat_ids = "compartmentalizer", "forgiving", "draconian"

	gen_skip = num_gens÷100

	fig, axs = plt.subplots(1, 2, figsize = (10, 5))
	ax = axs[1]
	[ax.plot(1:gen_skip:num_gens, strategy_freqs_array[1:gen_skip:num_gens,x], label=strat_ids[x]) for x in 1:3]
	ax.plot(1:gen_skip:num_gens, total_cooperation[1:gen_skip:num_gens], ls = "--", label="cooperation")
	ax.set_xlabel("time (Moran generations)")
	ax.set_ylabel("frequency")
	ax.set_ylim([0,1])
	ax.legend(loc=2)

	ax = axs[2]
	[ax.plot(1:gen_skip:num_gens, fitness_means_array[1:gen_skip:num_gens,x], label=strat_ids[x]) for x in 1:3]
	ax.plot(1:gen_skip:num_gens, fitness_means[1:gen_skip:num_gens], ls = "--", label="overall mean")
	ax.set_xlabel("time (Moran generations)")
	ax.set_ylabel("mean fitness by type")
	ax.legend(loc=2)

	plt.suptitle("random equal set membership, b = $b, c = $c, N=$N, M=$M, K=$K, w=$w, u_s=$u_s, u_p=$u_p, u_a=$u_p")
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	display(fig)
	plt.savefig("figures/simulation_run_redo_$(pci).pdf")
end
