# reputation_sets
Studying the evolution of reputation tracking mechanisms when individuals interact in multiple social contexts.

# basic usage notes
Install Julia from https://julialang.org/downloads/.
Ensure that the path to the module directory is in your Julia LOAD_PATH variable.
The easiest way to do this is to nano ~/.julia/config/startup.jl and add the line:
	push!(LOAD_PATH, "/path/to/reputation_sets/module")
Julia scripts can be run from the command line via julia src/name_of_script.jl
or by opening a Julia REPL and running the command:
	include("src/name_of_script.jl")
