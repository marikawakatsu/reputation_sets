#!/usr/bin/env julia

## parallel_lending.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Parallelized implementation of ReputationSets.

using Random, Statistics
using ReputationSets
using Distributed
using Revise
using ArgParse
using CSV
using Dates
using DataFrames
import JSON

function read_parameters(defpars::Dict{String,Any},
    inputfile = nothing)

    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey])
            if isprimitivetype(pars[parkey]["type"]) ||
                pars[parkey]["type"] == String
                T = pars[parkey]["type"]
            end
        else
            # default type is Float64
            T = Float64
        end
        #println(parkey, T)
        if T <: Int
            convertf = (val)->round(T, val)
        else
            convertf = (val)->convert(T, val)
        end

        # use defpars for list of usable parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r)->b.^r
                    pop!(valr, "log")
                else
                    rf = (r)->r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val)) || typeof(val) == String || typeof(val) == Bool
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

function main(args)

    s = ArgParseSettings(description =
        "run ReputationSets simulations across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type = Int64
            default = max(round(Int, Sys.CPU_THREADS / 2), 1)
        "--input"
            default = nothing
        #"--output"
        #    default=nothing
    end
    parsed_args = parse_args(args, s)

    defpars = Dict{String,Any}([
        "N"     => Dict("value" => 100, "type" => Int64),
        "M"     => Dict("value" => 3, "type" => Int64),
        "K"     => Dict("value" => 2, "type" => Int64),
		"b"     => Dict("value" => 1.0,     "type" => Float64),
		"c"     => Dict("value" => 0.1,     "type" => Float64),
        "δ"     => Dict("value" => 0.0,     "type" => Float64),
		"ϵ"     => Dict("value" => 0.1,     "type" => Float64),
		"w"     => Dict("value" => 0.01,     "type" => Float64),
		"u_s"     => Dict("value" => 0.01,     "type" => Float64),
		"u_p"     => Dict("value" => 0.01,     "type" => Float64),
		"u_a"     => Dict("value" => 0.01,     "type" => Float64),
        "num_samples" => Dict("value" => 5, "type" => Int64),
        "num_trials" => Dict("value" => 10, "type" => Int64),
        "output" => Dict("value" => "output/test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS / 2)))
    wpool = WorkerPool(workers())
    #extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)
    #@everywhere workers() push!(LOAD_PATH, $extradir)
    [@everywhere workers() push!(LOAD_PATH, $x) for x in extradir]
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using Statistics))
    @everywhere workers() eval(:(using ReputationSets))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))
    results = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))

    @everywhere function run_worker(inputs, results)
        # save trial number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.GLOBAL_RNG.seed))

        while true
            pard = take!(inputs)
            pard = merge(pard, seed)

            N = pard["N"]
			M = pard["M"]
			K = pard["K"]

			b = pard["b"]
			c = pard["c"]

			δ = pard["δ"]
			ϵ = pard["ϵ"]
			w = pard["w"]
			u_s = pard["u_s"]
			u_p = pard["u_p"]
			u_a = pard["u_a"]

			num_samples = pard["num_samples"]

            output = pard["output"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

			sets = equal_sets(N, M, K)
			game = Game(b, c, δ, ϵ, w, u_s, u_p, u_a, "db")
			pop = Population(sets, game, false)

			sampling_interval = N
			total_interactions = 2.0*sum([length(x) for x in sets.set_pairs])

			total_cooperation = Float64[]
			fitness_means = Float64[]
			strategy_freqs = Array{Float64, 1}[]
			strat_fitness_means = Array{Float64, 1}[]

			evolve!(pop, N^2)

			set_reputation_means = Float64[]
			reputation_means = Float64[]
			reputation_stds = Float64[]

			for g in 1:num_samples
				evolve!(pop, sampling_interval)
				push!(total_cooperation, sum(pop.prev_actions)/total_interactions)
				push!(fitness_means, mean(pop.fitnesses))
				tmp_reputations = [Int64[] for i in 1:pop.sets.N]
				tmp_set_reputations = [Int64[] for i in 1:pop.sets.M]
				for k in 1:pop.sets.M
					for j in pop.sets.set_members[k]
						if pop.sets.h[j,k]
							push!(tmp_reputations[j], pop.reputations[k, j])
							push!(tmp_set_reputations[k], pop.reputations[k, j])
						end
					end
				end
				set_rep_mean = [mean(x) for x in tmp_set_reputations]
				rep_mean = [mean(x) for x in tmp_reputations]
				if K > 1
					rep_std = mean([std(x) for x in tmp_reputations])
				else
					rep_std = 0.0
				end
				push!(reputation_means, mean(rep_mean))
				push!(reputation_stds, mean(rep_std))
				push!(set_reputation_means, mean(set_rep_mean))
				gen_freqs = zeros(Float64, 3)
				[gen_freqs[pop.strategies[i]+1] += 1.0/pop.sets.N for i in 1:N]
				push!(strategy_freqs, gen_freqs)
				gen_means = zeros(Float64, 3)
				[gen_means[x+1] += mean(pop.fitnesses[pop.strategies .== x]) for x in 0:2]
				push!(strat_fitness_means, gen_means)
			end

			strategy_freqs_array = zeros(Float64, num_samples, 3)
			fitness_means_array = zeros(Float64, num_samples, 3)
			for g in 1:num_samples
				strategy_freqs_array[g,:] = strategy_freqs[g]
				fitness_means_array[g,:] = strat_fitness_means[g]
			end
            #mean_freqs = [sum(indv==x for indv in pop.lattice)*1.0/N^2 for x in types_to_compete]
            pard["strat_fitness_means"] = fitness_means_array
			pard["fitness_means"] = fitness_means
			pard["strategy_freqs"] = strategy_freqs_array
			pard["total_cooperation"] = total_cooperation
			pard["set_reputation_means"] = set_reputation_means
			pard["reputation_means"] = reputation_means
			pard["reputation_stds"] = reputation_stds

            # return data to master process
            put!(results, pard)
        end
    end

    total_time_start = now()

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        #println(pard)
        println("--- queueing --- ")
        foreach(k->print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_trials"]
            nruns += 1
            rpard = copy(pard)
            rpard["rep"] = rep
            rpard["nrun"] = nruns
            put!(inputs, rpard)
        end
    end

    # start workers running on parameter sets in inputs
    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end

    # create output file name and data table
    output = pars["output"][1]
    println(output)
    file = occursin(r"\.csv$", output) ? output : output * ".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "strat_fitness_means", "fitness_means", "strategy_freqs", "total_cooperation", "seed1", "seed2", "seed3", "seed4"]...)
    dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    # grab results and output to CSV
    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        nrun = pop!(resd, "nrun")

        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end
    total_time_stop = now()
    total_time = Dates.canonicalize(Dates.CompoundPeriod(round(total_time_stop - total_time_start, Dates.Second(1))))
    println("total time elapsed: $total_time")
end

#main(ARGS)

main(["--input", "submit/simplified_M_K.json"])
