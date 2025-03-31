import os
from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, outdir_root, nvars, varlabels, varlabels_simple, rundir_name
from defaults import param_str, gamma, startidx_for_rom, basisname, params_train
from defaults import params_test, dt, error_lims
from defaults import FONTSIZE_LEGEND_LINE, FONTSIZE_TICKLABELS_LINE, FONTSIZE_TITLE_LINE, FONTSIZE_AXISLABEL_LINE
mpl.rc("axes", titlesize=FONTSIZE_TITLE_LINE)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_LINE)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_LINE)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.error_utils import calc_error_norms
from pschwarz.prom_utils import load_pod_basis, load_reduced_data


# ----- START USER INPUTS -----

domx = 2
domy = 2
overlap_list = [0, 10, 20]
schwarztype = "additive"
algos = ["LSPG"] * domx * domy

problem = "2d_swe"
varplot = 0
modelist = [40, 60, 80]

# problem = "2d_burgers"
# varplot = 0
# modelist = [100, 150, 200]

# problem = "2d_euler"
# varplot = 4
# modelist = [150, 200, 250]

plotcolors = ["m", "b", "g"]
plotstyles = ["-", "--", ":"]
legend_loc = "upper left"

mismatch_nan = True
gridlines = True

# ----- END USER INPUTS -----

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

basis_root = os.path.join(basedir, "pod_bases", basisname[problem], meshname, f"{domx}x{domy}") + "/overlap{overlap}/" + fluxorder

startlist = [startidx_for_rom[problem], 0]

for overlap_idx, overlap in enumerate(overlap_list):
    basisdir = basis_root.format(overlap=overlap)
    basis, center, norm = load_pod_basis(
        basisdir,
        return_basis=True,
        return_center=True,
        return_norm=True,
        nmodes=[np.max(modelist)]*domx*domy,
    )

    meshdir_decomp = os.path.join(basedir, "meshes", meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder, "full")

    for mode_idx, nmodes in enumerate(modelist):

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
                basedir, "runs", "rom_decomp",
                meshname, f"{domx}x{domy}",
                f"overlap{overlap}", fluxorder,
                rundir_name[problem].format(param=param),
                schwarztype, "_".join([f"{algo}{nmodes}" for algo in algos]),
            )

            rom_data = load_reduced_data(
                romdir,
                "state_snapshots",
                nvars[problem],
                meshdir_decomp,
                basisdir,
                "basis",
                "center",
                "norm",
                [nmodes]*domx*domy,
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

        if problem == "2d_burgers":
            artist, = ax.loglog(param_list, error_params, color=plotcolors[mode_idx], linestyle=plotstyles[overlap_idx])
        else:
            artist, = ax.semilogy(param_list, error_params, color=plotcolors[mode_idx], linestyle=plotstyles[overlap_idx])

        if overlap_idx == 0:
            artist_list.append(copy(artist))
            legend_labels.append(f"M = {nmodes}")

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


plt.tight_layout()
outfile = os.path.join(outdir, f"{varlabels_simple[problem][varplot]}_error_rom_decomp_spacetime.png")
print(f"Saving image to {outfile}")
plt.savefig(outfile)
plt.close(fig)

print("Finished")