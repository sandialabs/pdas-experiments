import os
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, rundir_name, param_str, params_train, params_test, outdir_root
from defaults import FONTSIZE_LEGEND_LINE, FONTSIZE_TICKLABELS_LINE, FONTSIZE_TITLE_LINE, FONTSIZE_AXISLABEL_LINE
mpl.rc("axes", titlesize=FONTSIZE_TITLE_LINE+2)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_LINE+2)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_LINE+2)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_LINE+2)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_LINE+2)

from pschwarz.data_utils import read_runtimes


# ----- START USER INPUTS ------

domx = 2
domy = 2
overlap_list = [0, 10, 20, 30]
schwarztype = "additive"

problem = "2d_swe"
nmodes = 80
plot_legend = True
legend_loc = "upper left"

# problem = "2d_burgers"
# nmodes = 150
# plot_legend = False

# problem = "2d_euler"
# nmodes = 200
# plot_legend = False

plotcolors = ["m", "b", "g", "c"]

gridlines = True
ylim = [1, 10]
yticks = [2, 4, 6, 8, 10]
yticklabels = [r"    2", r"    4", r"    6", r"    8", r"   10"]

# ----- END USER INPUTS -----

outdir = os.path.join(outdir_root, problem)
if not os.path.isdir(outdir):
    os.mkdir(outdir)

basedir = os.path.join(basedir_root, problem)

param_list = list(np.sort(params_test[problem] + params_train[problem]))
nparams = len(param_list)
noverlap = len(overlap_list)

fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
artist_list = []
legend_labels = []

runtime_arr = np.zeros((noverlap, nparams), dtype=np.float64)
for overlap_idx, overlap in enumerate(overlap_list):
    print(f"Overlap {overlap}")

    iters_arr = np.zeros(nparams, dtype=np.float64)

    for param_idx, param in enumerate(param_list):

        romdir = os.path.join(
            basedir, "runs", "rom_decomp",
            meshname, f"{domx}x{domy}",
            f"overlap{overlap}", fluxorder,
            rundir_name[problem].format(param=param), schwarztype,
            "_".join([f"LSPG{nmodes}"] * domx * domy),
        )

        runtime, iters, subiters = read_runtimes(romdir, "runtime")
        if (iters[0] == 0) or (subiters[0] == 0):
            iters_arr[param_idx] = np.nan
            runtime_arr[overlap_idx, param_idx] = np.nan
        else:
            iters_arr[param_idx] = subiters[0] / iters[0]
            runtime_arr[overlap_idx, param_idx] = runtime[0]

    if problem == "2d_burgers":
        ax.semilogx(param_list, iters_arr, color=plotcolors[overlap_idx])
    else:
        ax.plot(param_list, iters_arr, color=plotcolors[overlap_idx])

ax.set_ylim(ylim)
ax.set_ylabel("Average Schwarz iterations")
if plot_legend:
    ax.legend([r"$N_o$ = " + f"{overlap}" for overlap in overlap_list], loc=legend_loc)

if problem == "2d_burgers":
    ax.set_xticks(param_list)
    ax.tick_params(axis='x', which='minor', bottom=False)
else:
    ax.tick_params(labelbottom=False)

if problem == "2d_burgers":
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.tick_params(axis="y", which="major")
    ax.tick_params(axis="x", which="minor")
    ax.tick_params(labelbottom=False)
else:
    pass

ax.set_yticks(yticks, labels=yticklabels)

if gridlines:
    ax.grid(visible=True, which="major", axis="both", alpha=0.5, color="gray")

plt.tight_layout()
outfile = os.path.join(outdir, f"iterations_vs_overlap_k{nmodes}")
outfile += ".png"
print(f"Saving image to {outfile}")
plt.savefig(outfile)

plt.close(fig)

# ----- normalized runtime cost -----

fig, ax = plt.subplots(1, 1, figsize=(6.4, 5.2))

runtime_arr /= np.max(runtime_arr)
for overlap_idx, overlap in enumerate(overlap_list):
    if problem == "2d_burgers":
        ax.semilogx(param_list, runtime_arr[overlap_idx, :], color=plotcolors[overlap_idx])
    else:
        ax.plot(param_list, runtime_arr[overlap_idx, :], color=plotcolors[overlap_idx])

if problem == "2d_burgers":
    ax.set_xticks(param_list)
    ax.tick_params(axis='x', which='minor', bottom=False)
else:
    ax.set_xticks(param_list, minor=False)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

for param in params_test[problem]:
    ax.get_xticklabels()[param_list.index(param)].set_color("red")

if problem == "2d_burgers":
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.tick_params(axis="y", which="major")
    ax.tick_params(axis="x", which="minor")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
else:
    ax.tick_params(axis="both", which="major")

ax.set_xlabel(param_str[problem])
ax.set_ylabel("Normalized runtime")
ax.set_ylim([0, 1])
ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

if gridlines:
    ax.grid(visible=True, which="major", axis="both", alpha=0.5, color="gray")

plt.tight_layout()
outfile = os.path.join(outdir, f"normalized_cost_k{nmodes}")
outfile += ".png"
print(f"Saving image to {outfile}")
plt.savefig(outfile)

plt.close(fig)

print("Finished")
