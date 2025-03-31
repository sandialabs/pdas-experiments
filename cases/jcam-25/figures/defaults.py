import os

import matplotlib as mpl


# for contours
FONTSIZE_LEGEND_CONTOUR = 28
FONTSIZE_TITLE_CONTOUR = 32
FONTSIZE_AXISLABEL_CONTOUR = 28
FONTSIZE_TICKLABELS_CONTOUR = 26

# for line plots
FONTSIZE_LEGEND_LINE = 18
FONTSIZE_TITLE_LINE = 24
FONTSIZE_AXISLABEL_LINE = 20
FONTSIZE_TICKLABELS_LINE = 18

mpl.rc("font", family="serif")
mpl.rc("figure", facecolor="w", dpi=300)
mpl.rc("text", usetex=False)
mpl.rc("text.latex", preamble=r"\usepackage{amsmath}")

basedir_root = "<...>/pdas-experiments/cases/jcam-25"
outdir_root = "<...>/pdas-experiments/cases/jcam-25/figures/figs"

meshname = "300x300"
fluxorder = "firstorder"

if not os.path.isdir(outdir_root):
    os.mkdir(outdir_root)

bound_colors = ["r", "b", "c", "m"]

# for 2d_euler
gamma = 1.4

nvars = {
    "2d_swe": 3,
    "2d_burgers": 2,
    "2d_euler": 4,
}
varlabels = {
    "2d_swe": ["Water Height", "X-momentum", "Y-momentum"],
    "2d_burgers": ["X-velocity", "Y-velocity"],
    "2d_euler": ["Density", "X-momentum", "Y-momentum", "Energy", "Pressure"],
}
varlabels_simple = {
    "2d_swe": ["water_height", "xmomentum", "ymomentum"],
    "2d_burgers": ["xvel", "yvel"],
    "2d_euler": ["density", "xmomentum", "ymomentum", "energy", "pressure"],
}
nlevels = {
    "2d_swe": [97, 41, 41],
    "2d_burgers": [81, 81],
    "2d_euler": [29, 41, 41, 15, 61],
}
skiplevels = {
    "2d_swe": [32, 8, 8],
    "2d_burgers": [16, 16],
    "2d_euler": [2, 2, 2, 2, 12],
}
contourbounds = {
    "2d_swe": [[1.0, 1.024], [-0.05, 0.05], [-0.05, 0.05]],
    "2d_burgers": [[0.0, 0.5], [0.0, 0.5]],
    "2d_euler": [[0.1, 1.5], [-0.5, 0.5], [-0.5, 0.5], [0.25, 3.75], [0, 1.5]],
}
rundir_name = {
    "2d_swe": "coriolis{param}_gravity9.8_pulseMagnitude0.125_pulseX1.0_pulseY1.0",
    "2d_burgers": "diffusion{param}_pulseMagnitude0.5_pulseSpread0.075_pulseX-0.5_pulseY-0.4",
    "2d_euler": "gamma1.4_riemannTopRightPressure{param}_riemannTopRightXVel0.0_riemannTopRightYVel0.0_riemannTopRightDensity1.5_riemannBotLeftPressure0.029",
}
param_str = {
    "2d_swe": r"$\mu$",
    "2d_burgers": "D",
    "2d_euler": "p$_4$",
}
plotskip = {
    "2d_swe": 20,
    "2d_burgers": 6,
    "2d_euler": 5,
}
startidx_for_rom = {
    "2d_swe": 0,
    "2d_burgers": 0,
    "2d_euler": 10,
}
basisname = {
    "2d_swe": "coriolis_0p0_to_n4p0",
    "2d_burgers": "diffusion_0p0001_to_0p001_log",
    "2d_euler": "topRightPress_0p5_to_1p5",
}
params_train = {
    "2d_swe": [0.0, -1.0, -2.0, -3.0, -4.0],
    "2d_burgers": [0.0001, 0.000175, 0.000325, 0.00055, 0.001],
    "2d_euler": [0.5, 0.75, 1.0, 1.25, 1.5],
}
params_test = {
    "2d_swe": [-0.5, -1.5, -2.5, -3.5],
    "2d_burgers": [0.000135, 0.00025, 0.000425, 0.00075],
    "2d_euler": [0.625, 0.875, 1.125, 1.375],
}
dt = {
    "2d_swe": 0.01,
    "2d_burgers": 0.05,
    "2d_euler": 0.005,
}
error_lims = {
    "2d_swe": [1e-5, 2],
    "2d_burgers": [1e-5, 2],
    "2d_euler": [1e-5, 2],
}
tickvals = {
    "2d_swe": [-4.9833, -2.5, 0, 2.5, 4.9833],
    "2d_burgers": [-0.9967, -0.5, 0.0, 0.5, 0.9967],
    "2d_euler": [0.00167, 0.25, 0.5, 0.75, 0.9983],
}
ticklabels = {
    "2d_swe": ["-5", "-2.5", "0", "2.5", "5"],
    "2d_burgers": ["-1.0", "-0.5", "0.0", "0.5", "1.0"],
    "2d_euler": ["0.0", "0.25", "0.5", "0.75", "1.0"],
}