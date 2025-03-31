import os

from defaults import FONTSIZE_AXISLABEL_CONTOUR, FONTSIZE_TICKLABELS_CONTOUR, FONTSIZE_TITLE_CONTOUR, FONTSIZE_LEGEND_CONTOUR
from defaults import meshname, fluxorder
from defaults import basedir_root, outdir_root, nvars, nlevels, varlabels, skiplevels
from defaults import contourbounds, rundir_name, gamma, tickvals, ticklabels
import matplotlib as mpl
mpl.rc("axes", titlesize=FONTSIZE_TITLE_CONTOUR)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_CONTOUR)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_CONTOUR)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_CONTOUR)

from pschwarz.data_utils import load_unified_helper, euler_calc_pressure
from pschwarz.vis_utils import plot_contours


# ----- START USER INPUTS -----

problem = "2d_swe"
param = -2.0
varplot = 0
dt = 0.01
tm = 5.0
tf = 10.0

# problem = "2d_burgers"
# param = 0.000325
# varplot = 0
# dt = 0.05
# tm = 3.75
# tf = 7.0

# problem = "2d_euler"
# param = 1.0
# varplot = 4
# dt = 0.005
# tm = 0.45
# tf = 0.9

draw_colorbar = True

# ----- END USER INPUTS -----

outdir = os.path.join(outdir_root, problem, "contours", f"fom_{param}")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

basedir = os.path.join(basedir_root, problem)

meshdir = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")
datadir = os.path.join(basedir, "runs", "fom", meshname, fluxorder, rundir_name[problem].format(param=param))

meshlist, datalist = load_unified_helper(
    meshdirs=meshdir,
    datadirs=datadir,
    nvars=nvars[problem],
    dataroot="state_snapshots",
)

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

# plot first frame with colorbar
plot_contours(
    vardata,
    meshlist=meshlist,
    datalist=[datalist[0][:, :, [0], :]],
    nvars=nvars[problem],
    dataroot="state_snapshots",
    nlevels=nlevels[problem][varplot],
    skiplevels=skiplevels[problem][varplot],
    contourbounds=contourbounds[problem][varplot],
    plotskip=1,
    show_time=True,
    time_start=0.0,
    dt=dt,
    varlabel=varlabels[problem][varplot],
    figdim_base=[8, 9],
    vertical=False,
    draw_colorbar=True,
    tickvals=tickvals[problem],
    ticklabels=ticklabels[problem],
    savefigs=True,
    outdir=outdir,
    out_extra_str="_t0",
)

# draw middle and last frame without colorbar
nsnaps = datalist[0].shape[2]
plot_contours(
    vardata,
    meshlist=meshlist,
    datalist=[datalist[0][:, :, [int(nsnaps / 2)], :]],
    nvars=nvars[problem],
    dataroot="state_snapshots",
    nlevels=nlevels[problem][varplot],
    skiplevels=skiplevels[problem][varplot],
    contourbounds=contourbounds[problem][varplot],
    plotskip=1,
    show_time=True,
    time_start=tm,
    dt=dt,
    varlabel=varlabels[problem][varplot],
    figdim_base=[8, 8],
    vertical=False,
    draw_colorbar=False,
    tickvals=tickvals[problem],
    ticklabels=ticklabels[problem],
    savefigs=True,
    outdir=outdir,
    out_extra_str="_tm",
)

plot_contours(
    vardata,
    meshlist=meshlist,
    datalist=[datalist[0][:, :, [-1], :]],
    nvars=nvars[problem],
    dataroot="state_snapshots",
    nlevels=nlevels[problem][varplot],
    skiplevels=skiplevels[problem][varplot],
    contourbounds=contourbounds[problem][varplot],
    plotskip=1,
    show_time=True,
    time_start=tf,
    dt=dt,
    varlabel=varlabels[problem][varplot],
    figdim_base=[8, 8],
    vertical=False,
    draw_colorbar=False,
    tickvals=tickvals[problem],
    ticklabels=ticklabels[problem],
    savefigs=True,
    outdir=outdir,
    out_extra_str="_tf",
)