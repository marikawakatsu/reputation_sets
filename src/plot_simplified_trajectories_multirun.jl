#!/usr/bin/env julia

## plot_simplified_trajectories_multirun.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from ReputationSets simulations.
## Look at type frequency dynamics and cooperation
## as a function of c, M, and K.

using CSV, PyPlot, Statistics, ReputationSets

# load simulation output as a dataframe
runs = CSV.read("output/simplified_M_K_redo_bugfix_long.csv")

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Int64, String},Array{Float64, 3}}()
rep_means = Dict{Tuple{Int64, String},Array{Float64, 2}}()
rep_stds = Dict{Tuple{Int64, String},Array{Float64, 2}}()
set_rep_means = Dict{Tuple{Int64, String},Array{Float64, 3}}()
coop_freqs = Dict{Tuple{Int64, String},Array{Float64, 2}}()

N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
M_vals = sort(unique(runs[:M]))
runs[:permitted_strategies] = string.(runs[:permitted_strategies])
permitted_strategies = (sort(unique(runs[:permitted_strategies])))

num_samples = sort(unique(runs[:num_samples]))[1]
num_trials = sort(unique(runs[:num_trials]))[1]

param_combs = collect(Base.product(M_vals, permitted_strategies))

for (pi, param_comb) in enumerate(param_combs)
    M, permitted_strategies = param_comb
    println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, num_trials, num_samples, 5)
    coop_freqs[param_comb] = zeros(Float64, num_trials, num_samples)
	rep_means[param_comb] = zeros(Float64, num_trials, num_samples)
	rep_stds[param_comb] = zeros(Float64, num_trials, num_samples)
	set_rep_means[param_comb] = zeros(Float64, num_trials, num_samples, M)

    tmp_runs = runs[(runs[:M] .== M) .& (runs[:permitted_strategies] .== permitted_strategies), :]
    for (runi, run) in enumerate(eachrow(tmp_runs[:,:]))
		replicate = run[:rep]
        strat_rows = split.(run[:strategy_freqs][2:end-1], ";")
        coop_rows = split.(run[:total_cooperation][2:end-1], r", ")
		rep_means_rows = split.(run[:reputation_means][2:end-1], r", ")
		rep_stds_rows = split.(run[:reputation_stds][2:end-1], r", ")
		#set_rep_means_rows = split.(run[:set_reputation_means][2:end-1], ";")
        for (ri, row) in enumerate(strat_rows)
            #println("$row")
            if ri > 1
                row = row[2:end]
            end
            freqs = parse.(Float64,String.(split(row, " ")))
            #println("$freqs")
            type_freqs[param_comb][runi,ri,:] = freqs
        end
        for (ri, row) in enumerate(coop_rows)
            coop = parse(Float64, row)
            coop_freqs[param_comb][runi,ri] = coop
        end
		for (ri, row) in enumerate(rep_means_rows)
			means = parse(Float64, row)
			rep_means[param_comb][runi,ri] = means
		end
		for (ri, row) in enumerate(rep_stds_rows)
			means = parse(Float64, row)
			rep_stds[param_comb][runi,ri] = means
		end
		# if M > 1
		# 	for (ri, row) in enumerate(set_rep_means_rows)
	    #         #println("$row")
	    #         if ri > 1
	    #             row = row[2:end]
	    #         end
	    #         freqs = parse.(Float64,String.(split(row, " ")))
	    #         #println("$freqs")
	    #         set_rep_means[param_comb][ri,:] = freqs
	    #     end
		# end

    end
end

flattened_type_freqs = Dict{Tuple{Int64, String},Array{Float64, 2}}()
flattened_rep_means = Dict{Tuple{Int64, String},Array{Float64, 1}}()
flattened_rep_stds = Dict{Tuple{Int64, String},Array{Float64, 1}}()
flattened_set_rep_means = Dict{Tuple{Int64, String},Array{Float64, 1}}()
flattened_coop_freqs = Dict{Tuple{Int64, String},Array{Float64, 1}}()

for (pi, param_comb) in enumerate(param_combs)
    M, permitted_strategies = param_comb
    #println("$param_comb")
    flattened_type_freqs[param_comb] = zeros(Float64, num_trials*num_samples, 5)
    flattened_coop_freqs[param_comb] = zeros(Float64, num_trials*num_samples)
	flattened_rep_means[param_comb] = zeros(Float64, num_trials*num_samples)
	flattened_rep_stds[param_comb] = zeros(Float64, num_trials*num_samples)
	#flattened_set_rep_means[param_comb] = zeros(Float64, num_trials*num_samples, M)
	for ri in 1:num_trials
		flattened_type_freqs[param_comb][(ri-1)*num_samples+1:ri*num_samples,:] = type_freqs[param_comb][ri,:,:]
		flattened_coop_freqs[param_comb][(ri-1)*num_samples+1:ri*num_samples] = coop_freqs[param_comb][ri,:,:]
		flattened_rep_means[param_comb][(ri-1)*num_samples+1:ri*num_samples] = rep_means[param_comb][ri,:,:]
		flattened_rep_stds[param_comb][(ri-1)*num_samples+1:ri*num_samples] = rep_stds[param_comb][ri,:,:]
	end
end

param_combs = reshape(param_combs, length(param_combs))
strat_ids = "compartmentalizer", "forgiving", "draconian", "cooperator", "defector"

plotcolors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]


to_skip = 100
for (pci, param_comb) in enumerate(param_combs)
    fig = plt.figure()

	M, permitted_strategies = param_comb
	if permitted_strategies == "all"
		ps = collect(1:3)
	else
		ps = parse.(Int64, split(permitted_strategies, ""))
	end

	title_string = "M = $M"

    [plt.plot(collect(0:(100*to_skip):(num_samples*num_trials-1)*100), flattened_type_freqs[param_comb][1:to_skip:end,x], label=strat_ids[x], c = plotcolors[x]) for x in ps]
    plt.plot(collect(0:(100*to_skip):(num_samples*num_trials-1)*100), flattened_coop_freqs[param_comb][1:to_skip:end], ls = "--", c = plotcolors[4], label="cooperation")
	plt.plot(collect(0:(100*to_skip):(num_samples*num_trials-1)*100), flattened_rep_means[param_comb][1:to_skip:end], ls = "--", c = plotcolors[5], label=L"reputation \mu")
	plt.plot(collect(0:(100*to_skip):(num_samples*num_trials-1)*100), flattened_rep_stds[param_comb][1:to_skip:end], ls = "--", c = plotcolors[6], label=L"reputation \sigma")
    plt.xlabel("time (Moran generations)")
	[plt.vlines(x*num_samples*100, 0, 1, linestyles="--") for x in 1:num_trials-1]
	plt.ylabel("frequency")
	plt.ylim([0,1])
	plt.xlim([0,num_samples*num_trials*100])
	plt.legend(loc=2)
	plt.title(title_string)
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	display(fig)
	plt.savefig("figures/long_sim_M_$(M)_strategies_$(permitted_strategies)_redo.pdf")
end

# to_skip = 10
# for (pci, param_comb) in enumerate(param_combs)
# 	#if param_comb[2] == "all"
# 		for ri in 1:num_trials
# 		    fig = plt.figure()
#
# 			M, permitted_strategies = param_comb
# 			if permitted_strategies == "all"
# 				ps = collect(1:3)
# 			else
# 				ps = parse.(Int64, split(permitted_strategies, ""))
# 			end
#
# 			title_string = "M = $M, run $ri"
#
# 		    [plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), type_freqs[param_comb][ri,1:to_skip:end,x], label=strat_ids[x], c = plotcolors[x]) for x in ps]
# 		    plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), coop_freqs[param_comb][ri,1:to_skip:end], ls = "--", c = plotcolors[4], label="cooperation")
# 			plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), rep_means[param_comb][ri,1:to_skip:end], ls = "--", c = plotcolors[5], label=L"reputation \mu")
# 			plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), rep_stds[param_comb][ri,1:to_skip:end], ls = "--", c = plotcolors[6], label=L"reputation \sigma")
# 			# if M > 1
# 			# 	for i in 1:M
# 			# 		plt.plot(collect(0:(100*to_skip):(num_samples-1)*100), set_rep_means[param_comb][ri,1:to_skip:end], ls = "--", c = plotcolors[5], label="set $i reputation mu")
# 			# 	end
# 			# end
# 			plt.xlabel("time (Moran generations)")
# 			plt.ylabel("frequency")
# 			plt.ylim([0,1])
# 			plt.xlim([0,num_samples*100])
# 			plt.legend(loc=2)
# 			plt.title(title_string)
# 			fig.tight_layout(rect=[0, 0.03, 1, 0.96])
# 			display(fig)
# 			plt.savefig("figures/longer_sim_M_$(M)_permitted_strategies_$(permitted_strategies)_run_$(ri)_redo.pdf")
# 		end
# 	#end
# end
