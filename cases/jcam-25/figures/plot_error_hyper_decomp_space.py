import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, outdir_root, nvars, varlabels, rundir_name, bound_colors
from defaults import gamma, startidx_for_rom, basisname, params_train, params_test, dt, tickvals, ticklabels
from defaults import FONTSIZE_TICKLABELS_CONTOUR, FONTSIZE_TITLE_CONTOUR, FONTSIZE_AXISLABEL_CONTOUR
mpl.rc("axes", titlesize=FONTSIZE_TITLE_CONTOUR)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_CONTOUR)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.error_utils import calc_error_norms
from pschwarz.prom_utils import load_pod_basis, load_reduced_data
from pschwarz.vis_utils import plot_contours

# ----- START USER INPUTS -----

domx = 2
domy = 2
overlap = 0
schwarztype = "additive"
algos = ["LSPGHyper"] * domx * domy

problem = "2d_swe"
varplot = 0
param = -0.5
trial_modes = 80
greedy_modes = 0
gpod_modes = 0
sampperc = 0.005
domrate = 10
samp_algo = "random"
seed_qdeim = False
weigher = "identity"
error_log = True
error_lims = [-2, -1]
nlevels = 101
skiplevels = 100

# problem = "2d_burgers"
# varplot = 0
# param = 0.000135
# trial_modes = 150
# greedy_modes = 0
# gpod_modes = 0
# sampperc = 0.0375
# domrate = 2
# samp_algo = "random"
# seed_qdeim = True
# weigher = "identity"
# error_log = True
# error_lims = [-4, -2]
# nlevels = 81
# skiplevels = 40

# problem = "2d_euler"
# varplot = 4
# param = 1.375
# trial_modes = 200
# greedy_modes = 0
# gpod_modes = 0
# sampperc = 0.05
# domrate = 1
# samp_algo = "random"
# seed_qdeim = False
# weigher = "identity"
# error_log = True
# error_lims = [-2, 0]
# nlevels = 81
# skiplevels = 40

mismatch_nan = True

# ----- END USER INPUTS -----

outdir = os.path.join(outdir_root, problem, "error", "space_error")
if not os.path.isdir(outdir):
    os.mkdir(outdir)

basedir = os.path.join(basedir_root, problem)
meshdir_mono = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")

param_list = list(np.sort(params_test[problem] + params_train[problem]))
nparams = len(param_list)

fig, ax = plt.subplots(1, 1)
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
if samp_algo != "random":
    enddir += f"modes_{greedy_modes}_"
enddir += f"samp_{sampperc}"
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
    spacenorm=False,
    relative=False,
    mismatch_nan=mismatch_nan,
)

plot_contours(
    vardata,
    meshdirs=[meshdir_decomp],
    datalist=[np.expand_dims(errorlist[0], axis=2)],
    nvars=nvars[problem],
    dataroot="state_snapshots",
    draw_colorbar=False,
    contourlog=error_log,
    colormap="Reds",
    nlevels=nlevels,
    skiplevels=skiplevels,
    contourbounds=error_lims,
    plotbounds=True,
    bound_colors=bound_colors,
    tickvals=tickvals[problem],
    ticklabels=ticklabels[problem],
    varlabel=varlabels[problem][varplot],
    figdim_base=[8, 8],
    vertical=False,
    savefigs=True,
    outdir=outdir,
    colorbar_extra_str=" Absolute Error",
    out_extra_str="_space_hyper_decomp",
)
