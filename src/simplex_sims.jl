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

simplex_size = 10

simplex_freqs = Array{Float64, 1}[]
for i in 0:1:simplex_size
	for j in 0:1:(simplex_size - i)
		push!(simplex_freqs, [1.0/simplex_size*i, 1.0/simplex_size*j])
	end
end
#simplex_freqs = permutedims(hcat(simplex_freqs...))

b = 1.0 # benefit to cooperating
c = 0.8 # cost to cooperating

N = 200 # population size
M = 1 # number of sets
K = 1

δ = 0.0 # how much do aggregators favor their own set membership?
ϵ = 0.1 # cutoff for aggregators
w = 3.0/N # selection strength
u_s = 0.0/N # mutation rate between strategies
u_p = 3.0/N # error rate in choosing action
u_a = 3.0/N # error rate in assigning reputation

err = (1 - u_p)*(1 - u_a) - u_p*u_a

critical_freq = c/(b*(err - u_a))

permitted_strats = [1,5,4]

fixed_strats = []

results = []

for initial_freqs in simplex_freqs
	println("$initial_freqs")
	sets = equal_sets(N, M, K)
	game = Game(b, c, δ, ϵ, w, u_s, u_p, u_a, "db")
	pop = Population(sets, game, permitted_strats, false)

	pop_indices = floor.(Int64, N*initial_freqs)
	pop_indices[2] += pop_indices[1]

	pop.strategies[1:pop_indices[1]] .= permitted_strats[1]
	pop.strategies[pop_indices[1]+1:pop_indices[2]] .= permitted_strats[2]
	pop.strategies[pop_indices[2]+1:N] .= permitted_strats[3]

	total_interactions = 2.0*sum([length(x) for x in sets.set_pairs])

	strat_counts = zeros(Int64, 5)
	[strat_counts[pop.strategies[i]] += 1 for i in 1:N]
	while all(strat_counts .!= N)
		evolve!(pop)
		strat_counts = zeros(Int64, 5)
		[strat_counts[pop.strategies[i]] += 1 for i in 1:N]
	end
	fixed_strat = findfirst(strat_counts .== N)
	push!(fixed_strats, fixed_strat)

	result = vcat(initial_freqs, fixed_strat)
	push!(results, result)
end

function return_colors(strat::Int64)
	if strat < 4
		return 1
	else
		return strat - 2
	end
end

# legend_elements = [Line2D([0], [0], marker="o", label="discriminator",
# 	markerfacecolor=1),
# 	Line2D([0], [0], marker="o", label="defector",
# 	markerfacecolor=2),
# 	Line2D([0], [0], marker="o", label="cooperator",
# 	markerfacecolor=3)]

# results = permutedims(hcat(results...))
# fig = plt.figure()
# plt.scatter(results[:,1],results[:,2], c = [return_colors(floor(Int64,x)) for x in results[:,3]])
# plt.xlabel("discriminator frequency")
# plt.ylabel("defector frequency")
# plt.legend(loc=3)
# plt.tight_layout()
# display(fig)

disc_results = permutedims(hcat(filter(x->x[3] < 4, results)...))
coop_results = permutedims(hcat(filter(x->x[3] == 4, results)...))
def_results = permutedims(hcat(filter(x->x[3] == 5, results)...))

color = ["dodgerblue", "orange", "rebeccapurple"]

#results = permutedims(hcat(results...))
fig = plt.figure()
plt.scatter(disc_results[:,1],disc_results[:,2], c = color[1], label="Disc")
plt.scatter(coop_results[:,1],coop_results[:,2], c = color[2], label="AllC")
plt.scatter(def_results[:,1],def_results[:,2], c = color[3], label="AllD")

plt.plot([0,critical_freq], [0, 1 - critical_freq], ls = "--", c = "k", label="threshold")

plt.xlabel("discriminator frequency")
plt.ylabel("defector frequency")
plt.legend(loc=1)
plt.title("N = $N, (w, u_p, u_a) = ($w, $u_p, $u_a), b = $b, c = $c")
plt.tight_layout()


display(fig)
