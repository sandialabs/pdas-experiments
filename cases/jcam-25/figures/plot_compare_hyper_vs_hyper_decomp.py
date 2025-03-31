import os

import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, nvars, nlevels, varlabels, skiplevels
from defaults import contourbounds, rundir_name, gamma, plotskip, startidx_for_rom
from defaults import basisname, bound_colors, outdir_root, tickvals, ticklabels
from defaults import FONTSIZE_AXISLABEL_CONTOUR, FONTSIZE_TICKLABELS_CONTOUR, FONTSIZE_TITLE_CONTOUR, FONTSIZE_LEGEND_CONTOUR
mpl.rc("axes", titlesize=FONTSIZE_TITLE_CONTOUR)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_CONTOUR)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_CONTOUR)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.prom_utils import load_reduced_data
from pschwarz.vis_utils import plot_contours

# ----- START USER INPUTS -----

domx = 2
domy = 2
overlap = 0
schwarztype = "additive"

problem = "2d_swe"
param = -0.5
varplot = 0
time_start = 0.0
dt = 0.01
trial_modes_mono = 80
greedy_modes_mono = 0
gpod_modes_mono = 0
sampperc_mono = 0.005
trial_modes_decomp = 80
greedy_modes_decomp = 0
gpod_modes_decomp = 0
sampperc_decomp = 0.005
domrate = 10
samp_algo = "random"
seed_qdeim = False
weigher = "identity"

# problem = "2d_burgers"
# param = 0.000135
# varplot = 0
# time_start = 0.0
# dt = 0.05
# trial_modes_mono = 150
# greedy_modes_mono = 0
# gpod_modes_mono = 0
# sampperc_mono = 0.0375
# trial_modes_decomp = 150
# greedy_modes_decomp = 0
# gpod_modes_decomp = 0
# sampperc_decomp = 0.0375
# domrate = 2
# samp_algo = "random"
# seed_qdeim = True
# weigher = "identity"

# problem = "2d_euler"
# param = 1.375
# varplot = 4
# time_start = 0.05
# dt = 0.005
# trial_modes_mono = 200
# greedy_modes_mono = 0
# gpod_modes_mono = 0
# sampperc_mono = 0.05
# trial_modes_decomp = 200
# greedy_modes_decomp = 0
# gpod_modes_decomp = 0
# sampperc_decomp = 0.05
# domrate = 2
# samp_algo = "random"
# seed_qdeim = False
# weigher = "identity"

# ----- END USER INPUTS -----

outdir = os.path.join(outdir_root, problem, "contours", f"compare_hyper_vs_hyper_decomp_{param}")
if not os.path.isdir(outdir):
    os.mkdir(outdir)

basedir = os.path.join(basedir_root, problem)

meshdir_mono = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")
datadir = os.path.join(basedir, "runs", "fom", meshname, fluxorder, rundir_name[problem].format(param=param))

meshlist, datalist = load_unified_helper(
    meshdirs=meshdir_mono,
    datadirs=datadir,
    nvars=nvars[problem],
    dataroot="state_snapshots",
)

# load monolithic ROM
basis_dir_mono = os.path.join(basedir, "pod_bases", basisname[problem], meshname, "1x1", fluxorder)
romdir_mono = os.path.join(
    basedir, "runs",
    "hyper", meshname, fluxorder, rundir_name[problem].format(param=param),
    f"LSPG{trial_modes_mono}", samp_algo, weigher
)
enddir = ""
if weigher != "identity":
    enddir += f"modes_{gpod_modes_mono}_"
enddir += f"samp_{sampperc_mono}"
if seed_qdeim:
    enddir += "_qdeim"
romdir_mono = os.path.join(romdir_mono, enddir)
rom_data_mono = load_reduced_data(
    romdir_mono,
    "state_snapshots",
    nvars[problem],
    meshdir_mono,
    basis_dir_mono,
    "basis",
    "center",
    "norm",
    trial_modes_mono,
)
datalist.append(rom_data_mono.copy())

# load decomposed ROM
meshdir_decomp = os.path.join(basedir, "meshes", meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder, "full")
basis_dir_decomp = os.path.join(basedir, "pod_bases", basisname[problem], meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder)
romdir_decomp = os.path.join(
    basedir, "runs", "hyper_decomp",
    meshname, f"{domx}x{domy}",
    f"overlap{overlap}", fluxorder,
    rundir_name[problem].format(param=param),
    schwarztype, "_".join([f"LSPGHyper{trial_modes_decomp}"]*domx*domy),
    samp_algo, weigher
)
enddir = ""
if weigher != "identity":
    enddir += f"modes_{greedy_modes_decomp}_"
enddir += f"samp_{sampperc_decomp}"
if seed_qdeim:
    enddir += "_qdeim"
if domrate > 0:
    enddir += f"_dom{domrate}"
romdir_decomp = os.path.join(romdir_decomp, enddir)
rom_data_decomp = load_reduced_data(
    romdir_decomp,
    "state_snapshots",
    nvars[problem],
    meshdir_decomp,
    basis_dir_decomp,
    "basis",
    "center",
    "norm",
    [trial_modes_decomp]*domx*domy,
    merge_decomp=True,
)
datalist.append(rom_data_decomp.copy())

if (problem == "2d_euler") and (varplot == 4):
    datalist = euler_calc_pressure(
        gamma,
        meshlist=meshlist*len(datalist),
        datalist=datalist,
        nvars=nvars[problem],
        dataroot="state_snapshots",
    )
    vardata = 0
else:
    vardata = varplot

plotstart = [startidx_for_rom[problem], 0, 0]
plot_contours(
    vardata,
    meshdirs=[meshdir_mono, meshdir_mono, meshdir_decomp],
    datalist=datalist,
    nvars=nvars[problem],
    dataroot="state_snapshots",
    plotlabels=["Monolithic FOM", f"Monolithic HPROM", f"Schwarz HPROM"],
    nlevels=nlevels[problem][varplot],
    skiplevels=skiplevels[problem][varplot],
    contourbounds=contourbounds[problem][varplot],
    plotskip=plotskip[problem],
    plotstart=plotstart,
    varlabel=varlabels[problem][varplot],
    figdim_base=[8, 9],
    vertical=False,
    plotbounds=True,
    bound_colors=bound_colors,
    tickvals=tickvals[problem],
    ticklabels=ticklabels[problem],
    savefigs=True,
    outdir=outdir,
    show_time=True,
    time_start=time_start,
    dt=dt,
)
