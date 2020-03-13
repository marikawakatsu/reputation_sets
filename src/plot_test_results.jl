#!/usr/bin/env julia

## plt_test_results.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from ReputationSets simulations.
## Look at type frequency dynamics and cooperation
## as a function of c, M, and K.

using CSV, PyPlot, Statistics, NetworkLending

# load simulation output as a dataframe
runs = CSV.read("output/reputation_vary_u.csv")
save_fig = false

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Int64, Int64, Float64, Float64, Float64},Array{Float64, 1}}()
coop_freqs = Dict{Tuple{Int64, Int64, Float64, Float64, Float64},Float64}()

N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
M_vals = sort(unique(runs[:M]))
K_vals = sort(unique(runs[:K]))
c_vals = sort(unique(runs[:c]))
ua_vals = sort(unique(runs[:u_a]))
up_vals = sort(unique(runs[:u_p]))

param_combs = collect(Base.product(M_vals, K_vals, c_vals, ua_vals, up_vals))

for (pi, param_comb) in enumerate(param_combs)
    M, K, c, u_a, u_p = param_comb
    #println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, 3)
    coop_freqs[param_comb] = 0.0
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
            type_freqs[param_comb] += freqs/size(tmp_runs, 1)/run[:num_samples]
        end
        for (ri, row) in enumerate(coop_rows)
            coop = parse(Float64, row)
            coop_freqs[param_comb] += coop/size(tmp_runs, 1)/run[:num_samples]
        end
    end
end

param_combs = reshape(param_combs, length(param_combs))

ind = collect(1:length(param_combs))

freqbar = [[type_freqs[pc][i] for pc in param_combs] for i in 1:3]

fig = plt.figure()

#for pc in param_combs
plt.bar(ind, freqbar[1], label="compartmentalizer")
plt.bar(ind, freqbar[2], bottom=freqbar[1], label="forgiving")
plt.bar(ind, freqbar[3], bottom=freqbar[2]+freqbar[1], label="draconian")
plt.title("type frequencies")
plt.xticks(ind, param_combs, rotation=90)
plt.xlabel("M, K, c")
plt.ylabel("frequency")
plt.legend(loc=1)
plt.tight_layout()
display(fig)

if save_fig
    plt.savefig("figures/prelim_type_frequencies_redo.pdf")
end

fig = plt.figure()
plt.bar(ind, [coop_freqs[pc] for pc in param_combs])
plt.xticks(ind, param_combs, rotation=90)
plt.ylabel("frequency")
plt.title("cooperation fraction")
plt.ylim([0,1])
plt.xlabel("M, K, c")
plt.tight_layout()
display(fig)

if save_fig
    plt.savefig("figures/prelim_coop_frequencies_redo.pdf")
end

#
# zmin, zmax, rmin, rmax = z_vals[1], z_vals[end], r_vals[1], r_vals[end]
# scaling = (zmax-zmin)/(rmax-rmin)
#
# for (ki, k) in enumerate(k_vals)
#     fig, axs = plt.subplots(4, length(d_vals), figsize=(12,12),
#         sharey="col", sharex="row")
#     for (di, d) in enumerate(d_vals)
#         freqs_im = zeros(length(r_vals), length(z_vals), 4)
#         for (zi, z) in enumerate(z_vals)
#             for (ri, r) in enumerate(r_vals)
#                 freqs_im[ri, zi, :] += type_freqs[k, z, d, r]
#             end
#         end
#         for i in 1:4
#             ax = axs[i,di]
#             im = ax.imshow(freqs_im[:,:,i], origin = "lower",
#                 aspect=scaling,
#                 vmin=0, vmax=1,
#                 extent = [zmin, zmax, rmin, rmax])
#             if di == 1
#                 ax.set_ylabel("r")
#             end
#             if i == 1
#                 ax.set_title("d = $d")
#             end
#         end
#         #fig.colorbar(im)
#     #    cbar = ax.cax.colorbar(im)
#     #    cbar = grid.cbar_axes[0].colorbar(im)
#     end
#     fig.suptitle("payback, k = $k")
#     #plt.subplots_adjust(right=0.8)
#     #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     #plt.colorbar(im, cax=cbar_ax)
#     #fig.colorbar(im)
#     fig.tight_layout(rect=[0, 0.03, 1, 0.96])
#     #plt.subplots_adjust(top=0.85)
#     display(fig)
#
#
#     fig, axs = plt.subplots(2, length(d_vals), figsize=(9, 6),
#         sharey="col", sharex="row")
#     for (di, d) in enumerate(d_vals)
#         c_freq_im = zeros(length(r_vals), length(z_vals))
#         p_freq_im = zeros(length(r_vals), length(z_vals))
#         for (zi, z) in enumerate(z_vals)
#             for (ri, r) in enumerate(r_vals)
#                 c_freq_im[ri, zi] = c_freqs[k, z, d, r]
#                 p_freq_im[ri, zi] = p_freqs[k, z, d, r]
#             end
#         end
#         ax = axs[1,di]
#         im = ax.imshow(p_freq_im, origin = "lower",
#             aspect=scaling,
#             vmin=0, vmax=1,
#             extent = [zmin, zmax, rmin, rmax])
#         if di == 1
#             ax.set_ylabel("r")
#         end
#         ax.set_title("d = $d")
#
#         ax = axs[2,di]
#         im = ax.imshow(c_freq_im, origin = "lower",
#             aspect=scaling,
#             vmin=0, vmax=1,
#             extent = [zmin, zmax, rmin, rmax])
#         ax.set_xlabel("z")
#         if di == 1
#             ax.set_ylabel("r")
#         end
#         #fig.colorbar(im)
#     #    cbar = ax.cax.colorbar(im)
#     #    cbar = grid.cbar_axes[0].colorbar(im)
#     end
#     fig.suptitle("payback, k = $k")
#     #plt.subplots_adjust(right=0.8)
#     #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#     #plt.colorbar(im, cax=cbar_ax)
#     #fig.colorbar(im)
#     fig.tight_layout(rect=[0, 0.03, 1, 0.96])
#     #plt.subplots_adjust(top=0.85)
#     display(fig)
# end
