#!/usr/bin/env julia

using PyPlot

z_vals = collect(0.01:0.01:10)
R = 1
z2_law = 1.0 ./ (z_vals .^2)

full_expression = R ./ (R^2 .+ z_vals .^2)

fig  = plt.figure()
plt.plot(z_vals, z2_law, label="point source")
plt.plot(z_vals, full_expression, label="disc")

ax=plt.gca()
ax.set_xlim([0,5])
ax.set_ylim([0,5])
plt.vlines(R,0,5, linestyle = "--", label="R = 1")
plt.xlabel("distance (multiples of R)")
plt.ylabel("function value")
plt.legend(loc = 1)
plt.tight_layout()

display(fig)
plt.savefig("moon_landing.png")
