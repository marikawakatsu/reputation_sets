#!/usr/bin/env julia

## plt_test_results.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from ReputationSets simulations.
## Look at type frequency dynamics and cooperation
## as a function of c, M, and K.

using CSV, PyPlot, Statistics, ReputationSets

# load simulation output as a dataframe
runs = CSV.read("output/reputation_long_trajectory_M_K.csv")

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Int64, Int64, Float64, Float64, Float64},Array{Float64, 2}}()
coop_freqs = Dict{Tuple{Int64, Int64, Float64, Float64, Float64},Array{Float64, 1}}()

N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
M_vals = sort(unique(runs[:M]))
K_vals = sort(unique(runs[:K]))
c_vals = sort(unique(runs[:c]))
ua_vals = sort(unique(runs[:u_a]))
up_vals = sort(unique(runs[:u_p]))

num_samples = sort(unique(runs[:num_samples]))[1]

param_combs = collect(Base.product(M_vals, K_vals, c_vals, ua_vals, up_vals))

for (pi, param_comb) in enumerate(param_combs)
    M, K, c, u_a, u_p = param_comb
    #println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, num_samples, 3)
    coop_freqs[param_comb] = zeros(Float64, num_samples)
    tmp_runs = runs[(runs[:M] .== M) .& (runs[:K] .== K) .& (runs[:c] .== c) .& (runs[:u_a] .== u_a) .& (runs[:u_p] .== u_p), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
        strat_rows = split.(run[:strategy_freqs][2:end-1], ";")
        coop_rows = split.(run[:total_cooperation][2:end-1], r", ")
        for (ri, row) in enumerate(strat_rows)
            #println("$row")
            if ri > 1
                row = row[2:end]
            end
            freqs = parse.(Float64,String.(split(row, " ")))
            #println("$freqs")
            type_freqs[param_comb][ri,:] = freqs
        end
        for (ri, row) in enumerate(coop_rows)
            coop = parse(Float64, row)
            coop_freqs[param_comb][ri] = coop
        end
    end
end

param_combs = reshape(param_combs, length(param_combs))
strat_ids = "compartmentalizer", "forgiving", "draconian"

to_skip = 10

for (pci, param_comb) in enumerate(param_combs)
    fig = plt.figure()

	M, K, c, u_a, u_p = param_comb

	title_string = "M = $M, K = $K, c = $c, u_a = $u_a, u_p = $u_p"

    [plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), type_freqs[param_comb][1:to_skip:end,x], label=strat_ids[x]) for x in 1:3]
    plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), coop_freqs[param_comb][1:to_skip:end], ls = "--", label="cooperation")
    plt.xlabel("time (Moran generations)")
	plt.ylabel("frequency")
	plt.ylim([0,1])
	plt.xlim([0,num_samples*100])
	plt.legend(loc=2)
	plt.title(title_string)
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	display(fig)
	plt.savefig("figures/long_sim_$(pci)_M_$(M)_K_$(K).pdf")
end
