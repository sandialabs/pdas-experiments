import os
from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, outdir_root, nvars, varlabels, varlabels_simple, rundir_name, param_str
from defaults import gamma, startidx_for_rom, basisname, params_train, params_test, dt, error_lims
from defaults import FONTSIZE_LEGEND_LINE, FONTSIZE_TICKLABELS_LINE, FONTSIZE_TITLE_LINE, FONTSIZE_AXISLABEL_LINE
mpl.rc("axes", titlesize=FONTSIZE_TITLE_LINE+2)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_LINE+2)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_LINE+2)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_LINE+2)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_LINE)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.error_utils import calc_error_norms
from pschwarz.prom_utils import load_pod_basis, load_reduced_data

# ----- START USER INPUTS -----

domx = 2
domy = 2
overlap = 0
schwarztype = "additive"
algos = ["LSPGHyper"] * domx * domy

problem = "2d_swe"
varplot = 0
trial_modes = 80
greedy_modes = 0
gpod_modes = 0
sampperc_list = [0.005]
domrate_list = [5, 10, 30]
samp_algo = "random"
seed_qdeim = False
weigher = "identity"
legend_loc = "center right"

# problem = "2d_burgers"
# varplot = 0
# trial_modes = 150
# greedy_modes = 0
# gpod_modes = 0
# sampperc_list = [0.0375]
# domrate_list = [1, 2, 3]
# samp_algo = "random"
# seed_qdeim = True
# weigher = "identity"
# legend_loc = "lower right"

# problem = "2d_euler"
# varplot = 4
# trial_modes = 200
# greedy_modes = 0
# gpod_modes = 0
# sampperc_list = [0.05]
# domrate_list = [1, 2, 3]
# samp_algo = "random"
# seed_qdeim = False
# weigher = "identity"
# legend_loc = "lower right"

plotcolors = ["m", "b", "g"]
plotstyles = ["-", "--", ":"]

gridlines = True
mismatch_nan = True

# ----- END USER INPUTS -----

spec_format = False
if len(sampperc_list) == 1:
    spec_format = True

outdir = os.path.join(outdir_root, problem, "error")
if not os.path.isdir(outdir):
    os.mkdir(outdir)

basedir = os.path.join(basedir_root, problem)
meshdir_mono = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")

param_list = list(np.sort(params_test[problem] + params_train[problem]))
nparams = len(param_list)

fig, ax = plt.subplots(1, 1, figsize=(6.4, 5.2))
artist_list = []
legend_labels = []

basis_root = os.path.join(basedir, "pod_bases", basisname[problem], meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder)
basisdir = basis_root.format(overlap=overlap)
basis, center, norm = load_pod_basis(
    basisdir,
    return_basis=True,
    return_center=True,
    return_norm=True,
    nmodes=[trial_modes]*domx*domy,
)

meshdir_decomp = os.path.join(basedir, "meshes", meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder, "full")

startlist = [startidx_for_rom[problem], 0]

for domrate_idx, domrate in enumerate(domrate_list):
    for samp_perc_idx, samp_perc in enumerate(sampperc_list):

        error_params = []
        for param_idx, param in enumerate(param_list):

            print(f"Parameter: {param}")

            # load FOM data
            fomdir = os.path.join(basedir, "runs", "fom", meshname, fluxorder, rundir_name[problem].format(param=param))
            meshlist, datalist = load_unified_helper(
                meshdirs=[meshdir_mono],
                datadirs=[fomdir],
                nvars=nvars[problem],
                dataroot="state_snapshots",
            )

            # load decomp ROM
            romdir = os.path.join(
                basedir, "runs", "hyper_decomp",
                meshname, f"{domx}x{domy}",
                f"overlap{overlap}", fluxorder,
                rundir_name[problem].format(param=param),
                schwarztype, "_".join([f"{algo}{trial_modes}" for algo in algos]),
                samp_algo, weigher
            )
            enddir = ""
            if weigher != "identity":
                enddir += f"modes_{greedy_modes}_"
            enddir += f"samp_{samp_perc}"
            if seed_qdeim:
                enddir += "_qdeim"
            if domrate > 0:
                enddir += f"_dom{domrate}"
            romdir = os.path.join(romdir, enddir)

            rom_data = load_reduced_data(
                romdir,
                "state_snapshots",
                nvars[problem],
                meshdir_decomp,
                basisdir,
                "basis",
                "center",
                "norm",
                [trial_modes]*domx*domy,
                basis_in=basis,
                center_in=center,
                norm_in=norm,
                merge_decomp=True,
            )
            datalist.append(rom_data.copy())

            if (problem == "2d_euler") and (varplot == 4):
                datalist = euler_calc_pressure(
                    gamma,
                    meshlist=meshlist,
                    datalist=datalist,
                    nvars=nvars[problem],
                    dataroot="state_snapshots",
                )
                vardata = 0
            else:
                vardata = varplot


            # compute error
            errorlist, _ = calc_error_norms(
                meshlist=meshlist,
                datalist=datalist,
                nvars=nvars[problem],
                dataroot="state_snapshots",
                dtlist=dt[problem],
                startlist=startlist,
                samplist=[1, 1],
                timenorm=True,
                spacenorm=True,
                relative=True,
                mismatch_nan=mismatch_nan,
            )

            error_params.append(errorlist[0][vardata])

        if spec_format:
            plotcolor = plotcolors[domrate_idx]
            linestyle = "-"
            legend_labels.append(r"$N_b$ = " + str(domrate))
        else:
            plotcolor = plotcolors[samp_perc_idx]
            linestyle = plotstyles[domrate_idx]
            legend_labels.append(r"N$_s$" + f" = {samp_perc * 100}%, " + r"N$_b$ = " + str(domrate))

        if problem == "2d_burgers":
            artist, = ax.loglog(param_list, error_params, color=plotcolor, linestyle=linestyle)
        else:
            artist, = ax.semilogy(param_list, error_params, color=plotcolor, linestyle=linestyle)

        artist_list.append(copy(artist))


# space-time error
ax.set_ylim(error_lims[problem])
ax.set_title(f"{varlabels[problem][varplot]}")
ax.set_xlabel(param_str[problem])
ax.set_ylabel(r"Relative $\ell^2$ error")

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

if gridlines:
    ax.grid(visible=True, which="major", axis="both", alpha=0.5, color="gray")

if spec_format:
    ax.legend(artist_list, legend_labels, loc=legend_loc)
    plt.tight_layout()
else:
    # Put a legend below current axis
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.0875 * len(sampperc_list),
                    box.width, box.height * (1.0 - 0.0875 * len(sampperc_list))])
    ax.legend(artist_list, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.2),
            ncol=len(domrate_list), fontsize=10)

outfile = os.path.join(outdir, f"{varlabels_simple[problem][varplot]}_error_hyper_decomp_spacetime.png")
print(f"Saving image to {outfile}")
plt.savefig(outfile)
plt.close(fig)

print("Finished")