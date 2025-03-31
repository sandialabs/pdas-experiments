import os
from copy import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import FONTSIZE_AXISLABEL_LINE, FONTSIZE_TICKLABELS_LINE, FONTSIZE_TITLE_LINE, FONTSIZE_LEGEND_LINE
from defaults import basedir_root, meshname, fluxorder, outdir_root, nvars, varlabels, basisname
from defaults import params_test, params_train, dt, rundir_name, param_str, gamma
from defaults import error_lims, varlabels_simple
import matplotlib as mpl
mpl.rc("axes", titlesize=FONTSIZE_TITLE_LINE)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_LINE)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_LINE)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.prom_utils import load_pod_basis, calc_projection
from pschwarz.error_utils import calc_error_norms

# ----- START USER INPUTS -----

domx = 2
domy = 2
overlap = 0
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

gridlines = True

# ----- END USER INPUTS -----

outdir = os.path.join(outdir_root, problem, "error")
if not os.path.isdir(outdir):
    os.mkdir(outdir)

basedir = os.path.join(basedir_root, problem)
meshdir_mono = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")
meshdir_decomp = os.path.join(basedir, "meshes", meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder, "full")

param_list = list(np.sort(params_test[problem] + params_train[problem]))
nparams = len(param_list)

# load monolithic basis
basisdir_mono = os.path.join(basedir, "pod_bases", basisname[problem], meshname, "1x1", fluxorder)
basis_mono, center_mono, norm_mono = load_pod_basis(
    basisdir_mono,
    return_basis=True,
    return_center=True,
    return_norm=True,
    nmodes=np.amax(modelist),
)

# load decomposed basis
basisdir_decomp = os.path.join(basedir, "pod_bases", basisname[problem], meshname, f"{domx}x{domy}", f"overlap{overlap}",  fluxorder)
basis_decomp, center_decomp, norm_decomp = load_pod_basis(
    basisdir_decomp,
    return_basis=True,
    return_center=True,
    return_norm=True,
    nmodes=[np.amax(modelist)]*domx*domy,
)

nparams = len(param_list)
nmodes = len(modelist)
ndomains = len(basis_decomp)
error_arr_scalar_mono   = np.zeros((nmodes, nparams), dtype=np.float64)
error_arr_scalar_decomp = np.zeros((nmodes, nparams), dtype=np.float64)
error_arr_time_mono  = None
error_arr_time_decomp = None

# loop cases
params_time_counter = 0
for param_idx, param in enumerate(param_list):

    print(f"Parameter: {param}")

    # load monolithic FOM data
    fomdir = os.path.join(basedir, "runs", "fom", meshname, fluxorder, rundir_name[problem].format(param=param))
    meshlist_mono, datalist_mono = load_unified_helper(
        meshdirs=meshdir_mono,
        datadirs=fomdir,
        nvars=nvars[problem],
        dataroot="state_snapshots",
    )

    meshlist = meshlist_mono * 3

    # loop modes
    for mode_idx, modes in enumerate(modelist):

        print(f"Modes: {modes}")

        datalist_mono_proj = calc_projection(
            modes,
            meshlist=meshlist_mono,
            datalist=datalist_mono,
            nvars=nvars[problem],
            dataroot="state_snapshots",
            basis_in=basis_mono,
            center=center_mono,
            norm=norm_mono,
        )


        datalist_decomp_proj = calc_projection(
            modes,
            meshlist=meshlist_mono,
            datalist=datalist_mono,
            nvars=nvars[problem],
            dataroot="state_snapshots",
            basis_in=basis_decomp,
            center=center_decomp,
            norm=norm_decomp,
            meshdir_decomp=meshdir_decomp,
            merge_decomp=True,
        )

        err_mesh = meshdir_mono


        if (problem == "2d_euler") and (varplot == 4):
            datalist = euler_calc_pressure(
                gamma,
                meshlist=meshlist,
                datalist=[datalist_mono[0], datalist_mono_proj[0], datalist_decomp_proj[0]],
                nvars=nvars[problem],
                dataroot="state_snapshots",
            )
            vardata = 0
        else:
            datalist = [datalist_mono[0], datalist_mono_proj[0], datalist_decomp_proj[0]]
            vardata = varplot

        errorlist_scalar, _ = calc_error_norms(
            meshlist=meshlist,
            datalist=datalist,
            nvars=nvars[problem],
            dataroot="state_snapshots",
            dtlist=dt[problem],
            samplist=[1, 1, 1],
            timenorm=True,
            spacenorm=True,
            relative=True,
            merge_decomp=True,
        )

        # store scalar errors
        error_arr_scalar_mono[mode_idx, param_idx] = errorlist_scalar[0][vardata].copy()
        error_arr_scalar_decomp[mode_idx, param_idx] = errorlist_scalar[1][vardata].copy()

# plot scalar
fig, ax = plt.subplots(figsize=(6.4, 5.2))

artists = []
for mode_idx, modes in enumerate(modelist):

    if problem == "2d_burgers":
        artist_mono, = ax.loglog(
            param_list, error_arr_scalar_mono[mode_idx, :],
            color=plotcolors[mode_idx],
            linestyle="-",
        )
        artist_decomp, = ax.loglog(
            param_list, error_arr_scalar_decomp[mode_idx, :],
            color=plotcolors[mode_idx],
            linestyle="--",
        )
    else:
        artist_mono, = ax.semilogy(
            param_list, error_arr_scalar_mono[mode_idx, :],
            color=plotcolors[mode_idx],
            linestyle="-",
        )
        artist_decomp, = ax.semilogy(
            param_list, error_arr_scalar_decomp[mode_idx, :],
            color=plotcolors[mode_idx],
            linestyle="--",
        )

    artists.append(copy(artist_mono))

ax.set_ylim(error_lims[problem])
ax.set_xlabel(param_str[problem])
ax.set_ylabel("Relative projection error")
ax.set_title(varlabels[problem][varplot])

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

legend_labels = [f"M = {modes}" for modes in modelist]
ax.legend(artists, legend_labels, loc="upper left", framealpha=1.0)

plt.tight_layout()
outfile = os.path.join(outdir, f"{varlabels_simple[problem][varplot]}_error_rom_decomp_proj_spacetime.png")
print(f"Saving image to {outfile}")
plt.savefig(outfile)
plt.close(fig)

print("Finished")
