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

using Revise, StatsBase, Statistics, PyPlot, Profile

using ReputationSets

b = 1.0 # benefit to cooperating
c = 0.4 # cost to cooperating

N = 50 # population size
M = 1 # number of sets
K = 1

δ = 0.0 # how much do aggregators favor their own set membership?
ϵ = 0.1 # cutoff for aggregators
w = 1.0/N # selection strength
u_s = 1.0/N # mutation rate between strategies
u_p = 3.0/N # error rate in choosing action
u_a = 3.0/N # error rate in assigning reputation

err = (1 - u_p)*(1 - u_a) - u_p*u_a

critical_freq = c/(b*(err - u_a))

permitted_strats = [0,3,4]

sets = equal_sets(N, M, K)
game = Game(b, c, δ, ϵ, w, u_s, u_p, u_a, "db")
pop = Population(sets, game, permitted_strats, false)

num_gens = 100000
gen_skip = num_gens÷1000


total_interactions = 2.0*sum([length(x) for x in sets.set_pairs])

total_cooperation = Float64[]
fitness_means = Float64[]
strategy_freqs = Array{Float64, 1}[]
strat_fitness_means = Array{Float64, 1}[]

evolve!(pop)

reputation_means = Float64[]
reputation_stds = Float64[]

for g in 1:num_gens
	#println("$g")
	evolve!(pop, 1)
	push!(total_cooperation, sum(pop.prev_actions)/total_interactions)
	push!(fitness_means, mean(pop.fitnesses))
	#tmp_reputations = zeros(pop.sets.N, pop.sets.K)
	tmp_reputations = [Int64[] for i in 1:pop.sets.N]
	for k in 1:pop.sets.M
		for j in pop.sets.set_members[k]
			if pop.sets.h[j,k]
				push!(tmp_reputations[j], pop.reputations[k, j])
			end
		end
	end
	rep_means = [mean(x) for x in tmp_reputations]
	rep_stds = [std(x) for x in tmp_reputations]
	gen_freqs = zeros(Float64, 5)
	[gen_freqs[pop.strategies[i]+1] += 1.0/pop.sets.N for i in 1:N]
	push!(strategy_freqs, gen_freqs)
	gen_means = zeros(Float64, 5)
	[gen_means[x+1] += mean(pop.fitnesses[pop.strategies .== x]) for x in 0:4]
	push!(strat_fitness_means, gen_means)
	push!(reputation_means, mean(rep_means))
	push!(reputation_stds, mean(rep_stds))
end

strategy_freqs_array = zeros(Float64, num_gens, 5)
fitness_means_array = zeros(Float64, num_gens, 5)
for g in 1:num_gens
	strategy_freqs_array[g,:] = strategy_freqs[g]
	fitness_means_array[g,:] = strat_fitness_means[g]
end

strat_ids = "compartmentalizer", "forgiving", "draconian", "cooperator", "defector"

fig, axs = plt.subplots(1, 1, figsize = (10, 5))
#ax = axs[1]
ax=plt.gca()
[ax.plot(1:gen_skip:num_gens, strategy_freqs_array[1:gen_skip:num_gens,x], label=strat_ids[x]) for x in permitted_strats .+ 1]
#ax.plot(1:gen_skip:num_gens, total_cooperation[1:gen_skip:num_gens], ls = "--", label="cooperation")
ax.set_xlabel("time (Moran generations)")
ax.set_ylabel("frequency")
ax.set_ylim([0,1])
ax.set_xlim([0,num_gens])
ax.hlines(critical_freq, 0, num_gens, linestyle = "--", label="critical frequency")
ax.legend(loc=2)

# ax = axs[2]
# [ax.plot(1:gen_skip:num_gens, fitness_means_array[1:gen_skip:num_gens,x], label=strat_ids[x]) for x in permitted_strats .+ 1]
# ax.plot(1:gen_skip:num_gens, fitness_means[1:gen_skip:num_gens], ls = "--", label="overall mean")
# ax.set_xlabel("time (Moran generations)")
# ax.set_ylabel("mean fitness by type")
# ax.legend(loc=2)
# #ax.savefig("figures/test_fitnesses_5.pdf")

plt.suptitle("random equal set membership, b = $b, c = $c, N=$N, M=$M, K=$K, w=$w, u_s=$u_s, u_p=$u_p, u_a=$u_p")
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
display(fig)
