import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from defaults import basedir_root, meshname, fluxorder, outdir_root, nvars, rundir_name, gamma, params_test
from defaults import startidx_for_rom, basisname, dt, error_lims
from defaults import FONTSIZE_LEGEND_LINE, FONTSIZE_TICKLABELS_LINE, FONTSIZE_TITLE_LINE, FONTSIZE_AXISLABEL_LINE
mpl.rc("axes", titlesize=FONTSIZE_TITLE_LINE)
mpl.rc("axes", labelsize=FONTSIZE_AXISLABEL_LINE)
mpl.rc("xtick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("ytick", labelsize=FONTSIZE_TICKLABELS_LINE)
mpl.rc("legend", fontsize=FONTSIZE_LEGEND_LINE-4)

from pschwarz.data_utils import read_runtimes, load_field_data, euler_calc_pressure
from pschwarz.prom_utils import load_reduced_data, load_pod_basis
from pschwarz.error_utils import calc_error_norms

# ----- START USER INPUTS -----

domx = 2
domy = 2
schwarztype = "additive"

problem = "2d_swe"
varplot = 0
param_list = params_test[problem]
overlap = 0
nmodes_list = [20, 40, 60, 80]
nmodes_list_decomp = [20, 40, 60, 80]
nmodes_hyper = 80
nmodes_hyper_decomp = 80
samppercs_hyper = [0.0025, 0.00375, 0.005, 0.00875, 0.01, 0.015, 0.025, 0.0375, 0.05]
samppercs_hyper_decomp = [0.005, 0.00875, 0.01, 0.015, 0.025, 0.0375, 0.05]
domrate = 10
sampalgo = "random"
seed_qdeim = False
weigher = "identity"
plot_legend = True

# problem = "2d_burgers"
# varplot = 0
# param_list = params_test[problem]
# overlap = 0
# nmodes_list = [100, 150, 200, 250, 300]
# nmodes_list_decomp = [100, 150, 200, 250, 300]
# nmodes_hyper = 150
# nmodes_hyper_decomp = 150
# samppercs_hyper = [0.01, 0.015, 0.025, 0.0375, 0.05, 0.0875, 0.1]
# samppercs_hyper_decomp = [0.01, 0.015, 0.025, 0.0375, 0.05, 0.0875, 0.1]
# domrate = 2
# sampalgo = "random"
# seed_qdeim = True
# weigher = "identity"
# plot_legend = False

# problem = "2d_euler"
# varplot = 4
# param_list = params_test[problem]
# overlap = 0
# nmodes_list = [200, 250, 300]
# nmodes_list_decomp = [150, 200, 250, 300]
# nmodes_hyper = 200
# nmodes_hyper_decomp = 200
# samppercs_hyper = [0.015, 0.025, 0.0375, 0.05, 0.0875, 0.1]
# samppercs_hyper_decomp = [0.015, 0.025, 0.0375, 0.05, 0.0875, 0.1]
# domrate = 2
# sampalgo = "random"
# seed_qdeim = False
# weigher = "identity"
# plot_legend = False

plotcolors = ["b", "r", "c", "m"]
legend_loc = "best"

xlim = [0.1, 1000]

mismatch_nan = True

plot_hyper = True
plot_hyper_decomp = True

gridlines = True

# ----- END USER INPUTS -----

algos = ["LSPG"] * domx * domy

outdir = os.path.join(outdir_root, problem)
if not os.path.isdir(outdir):
    os.mkdir(outdir)

if (problem == "2d_euler") and (varplot == 4):
    vardata = 0
else:
    vardata = varplot

def remove_exploded(err_in, speed_in, speed_decomp_in):
    err_out = []
    speed_out = []
    speed_decomp_out = []
    for val_idx, err_val in enumerate(err_in):
        if err_val < 0.1:
            err_out.append(err_val)
            speed_out.append(speed_in[val_idx])
            speed_decomp_out.append(speed_decomp_in[val_idx])
    return err_out, speed_out, speed_decomp_out

basedir = os.path.join(basedir_root, problem)
meshdir_mono = os.path.join(basedir, "meshes", meshname, "1x1", fluxorder, "full")
meshdir_decomp = os.path.join(basedir, "meshes", meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder, "full")

# load monolithic basis
basis_dir_mono = os.path.join(basedir, "pod_bases", basisname[problem], meshname, "1x1", fluxorder)
basis_mono, center_mono, norm_mono = load_pod_basis(
    basis_dir_mono,
    return_basis=True,
    return_center=True,
    return_norm=True,
    nmodes=max(nmodes_hyper, np.amax(nmodes_list)),
)
basis_dir_decomp = os.path.join(basedir, "pod_bases", basisname[problem], meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder)
basis_decomp, center_decomp, norm_decomp = load_pod_basis(
    basis_dir_decomp,
    return_basis=True,
    return_center=True,
    return_norm=True,
    nmodes=[max(nmodes_hyper, np.amax(nmodes_list_decomp))] * domx * domy,
)

fig, ax = plt.subplots(1, 1)
ax2 = ax.twiny()

for param_idx, param in enumerate(param_list):
    # FOM for relative runtime / error
    fomdir_mono = os.path.join(
        basedir, "runs", "fom",
        meshname, fluxorder,
        rundir_name[problem].format(param=param)
    )
    fomdata, _ = load_field_data(
        fomdir_mono,
        "state_snapshots",
        nvars[problem],
        meshdir=meshdir_mono,
        merge_decomp=False,
    )
    if (problem == "2d_euler") and (varplot == 4):
        fomdata_list = euler_calc_pressure(
            gamma,
            meshdirs=[meshdir_mono],
            datalist=[fomdata],
            nvars=nvars[problem],
            dataroot="state_snapshots",
        )
    else:
        fomdata_list = [fomdata]
    fom_time, _, _ = read_runtimes(fomdir_mono, "runtime")
    fom_time = fom_time[0]

    fomdir_decomp = os.path.join(
        basedir, "runs", "fom_decomp",
        meshname, f"{domx}x{domy}", f"overlap{overlap}", fluxorder,
        rundir_name[problem].format(param=param),
        schwarztype, "_".join(["FOM"]*domx*domy)
    )
    fom_time_decomp, _, _ = read_runtimes(fomdir_decomp, "runtime")
    fom_time_decomp = fom_time_decomp[0]

    # load monolithic PROM, various mode numbers
    romdata_list = []
    romdirs = []
    for mode_idx, nmodes in enumerate(nmodes_list):
        romdir = os.path.join(
            basedir, "runs", "rom",
            meshname, fluxorder,
            rundir_name[problem].format(param=param),
            f"LSPG{nmodes}",
        )
        time, _, _ = read_runtimes(romdir, "runtime")
        if time[0] == 0.0:
            continue

        romdata = load_reduced_data(
            romdir,
            "state_snapshots",
            nvars[problem],
            meshdir_mono,
            basis_dir_mono,
            "basis",
            "center",
            "norm",
            nmodes,
            basis_in=basis_mono,
            center_in=center_mono,
            norm_in=norm_mono,
        )

        romdirs.append(romdir)
        romdata_list.append(romdata.copy())

    nrom = len(romdirs)
    if (problem == "2d_euler") and (varplot == 4):
        romdata_list = euler_calc_pressure(
            gamma,
            meshdirs=[meshdir_mono] * nrom,
            datalist=romdata_list,
            nvars=nvars[problem],
            dataroot="state_snapshots",
        )

    # load decomposed PROM, various mode numbers
    romdecompdirs = []
    romdecompdata_list = []
    for mode_idx, nmodes in enumerate(nmodes_list_decomp):
        romdir = os.path.join(
            basedir, "runs", "rom_decomp",
            meshname, f"{domx}x{domy}",
            f"overlap{overlap}", fluxorder,
            rundir_name[problem].format(param=param),
            schwarztype,
            "_".join([f"LSPG{nmodes}"] * domx * domy),
        )
        time, _, _ = read_runtimes(romdir, "runtime")
        if time[0] == 0.0:
            continue

        romdata = load_reduced_data(
            romdir,
            "state_snapshots",
            nvars[problem],
            meshdir_decomp,
            basis_dir_decomp,
            "basis",
            "center",
            "norm",
            [nmodes] * domx * domy,
            basis_in=basis_decomp,
            center_in=center_decomp,
            norm_in=norm_decomp,
            merge_decomp=True,
        )

        romdecompdirs.append(romdir)
        romdecompdata_list.append(romdata.copy())

    nromdecomp = len(romdecompdata_list)
    if (problem == "2d_euler") and (varplot == 4):
        romdecompdata_list = euler_calc_pressure(
            gamma,
            meshdirs=[meshdir_mono] * nromdecomp,
            datalist=romdecompdata_list,
            nvars=nvars[problem],
            dataroot="state_snapshots",
        )

    # load monolithic HPROM, various sampling rates
    if plot_hyper:
        hyperdata_list = []
        hyperdirs = []
        for samp_idx, sampperc in enumerate(samppercs_hyper):
            hyper_dir = os.path.join(
                basedir, "runs", "hyper",
                meshname, fluxorder,
                rundir_name[problem].format(param=param),
                f"LSPG{nmodes_hyper}",
                sampalgo, weigher,
            )
            enddir = ""
            if weigher != "identity":
                enddir += f"modes_{nmodes_hyper}_"
            enddir += f"samp_{sampperc}"
            if seed_qdeim:
                enddir += "_qdeim"
            hyper_dir = os.path.join(hyper_dir, enddir)

            time, _, _ = read_runtimes(hyper_dir, "runtime")
            if time[0] == 0.0:
                continue

            hyperdata = load_reduced_data(
                hyper_dir,
                "state_snapshots",
                nvars[problem],
                meshdir_mono,
                basis_dir_mono,
                "basis",
                "center",
                "norm",
                nmodes_hyper,
                basis_in=basis_mono,
                center_in=center_mono,
                norm_in=norm_mono,
            )

            hyperdirs.append(hyper_dir)
            hyperdata_list.append(hyperdata.copy())

        nhyper = len(hyperdata_list)
        if (problem == "2d_euler") and (varplot == 4):
            hyperdata_list = euler_calc_pressure(
                gamma,
                meshdirs=[meshdir_mono] * nhyper,
                datalist=hyperdata_list,
                nvars=nvars[problem],
                dataroot="state_snapshots",
            )

    if plot_hyper_decomp:
        # load decomposed HPROM, various sampling rates
        hyperdecompdirs = []
        hyperdecompdata_list = []
        for samp_idx, sampperc in enumerate(samppercs_hyper_decomp):

            hyper_decomp_dir = os.path.join(
                basedir, "runs", "hyper_decomp",
                meshname, f"{domx}x{domy}",
                f"overlap{overlap}", fluxorder,
                rundir_name[problem].format(param=param),
                schwarztype,
                "_".join([f"LSPGHyper{nmodes_hyper_decomp}"] * domx * domy),
                sampalgo, weigher,
            )
            enddir = ""
            if weigher != "identity":
                enddir += f"modes_{nmodes_hyper_decomp}_"
            enddir += f"samp_{sampperc}"
            if seed_qdeim:
                enddir += "_qdeim"
            if domrate > 0:
                enddir += f"_dom{domrate}"
            hyper_decomp_dir = os.path.join(hyper_decomp_dir, enddir)

            if not os.path.isdir(hyper_decomp_dir):
                continue
            else:
                time, _, _ = read_runtimes(hyper_decomp_dir, "runtime")
                if time[0] == 0.0:
                    continue

            hyperdata = load_reduced_data(
                hyper_decomp_dir,
                "state_snapshots",
                nvars[problem],
                meshdir_decomp,
                basis_dir_decomp,
                "basis",
                "center",
                "norm",
                [nmodes_hyper_decomp] * domx * domy,
                basis_in=basis_decomp,
                center_in=center_decomp,
                norm_in=norm_decomp,
                merge_decomp=True,
            )

            hyperdecompdirs.append(hyper_decomp_dir)
            hyperdecompdata_list.append(hyperdata.copy())

        nhyperdecomp = len(hyperdecompdirs)
        if (problem == "2d_euler") and (varplot == 4):
            hyperdecompdata_list = euler_calc_pressure(
                gamma,
                meshdirs=[meshdir_mono] * nhyperdecomp,
                datalist=hyperdecompdata_list,
                nvars=nvars[problem],
                dataroot="state_snapshots",
            )

    # calculate errors and runtimes
    rom_errors, _ = calc_error_norms(
        meshdirs=[meshdir_mono] * (nrom + 1),
        datalist=fomdata_list + romdata_list,
        nvars=nvars[problem],
        dataroot="state_snapshots",
        dtlist=dt[problem],
        samplist=[1] * (nrom + 1),
        startlist=[startidx_for_rom[problem]] + [0] * nrom,
        timenorm=True,
        spacenorm=True,
        relative=True,
        mismatch_nan=mismatch_nan,
    )
    rom_errors_plot = [err[vardata] for err in rom_errors]
    rom_times, _, _ = read_runtimes(romdirs, "runtime")
    rom_speedup = [fom_time / time for time in rom_times]
    rom_speedup_v_decomp = [fom_time_decomp * domx * domy / time for time in rom_times]
    rom_errors_plot, rom_speedup, rom_speedup_v_decomp = remove_exploded(
        rom_errors_plot,
        rom_speedup,
        rom_speedup_v_decomp
    )

    rom_decomp_errors, _ = calc_error_norms(
        meshdirs=[meshdir_mono] * (1 + nromdecomp),
        datalist=fomdata_list + romdecompdata_list,
        nvars=nvars[problem],
        dataroot="state_snapshots",
        dtlist=dt[problem],
        samplist=[1] * (nromdecomp + 1),
        startlist=[startidx_for_rom[problem]] + [0] * nromdecomp,
        timenorm=True,
        spacenorm=True,
        relative=True,
        mismatch_nan=mismatch_nan,
    )
    rom_decomp_errors_plot = [err[vardata] for err in rom_decomp_errors]
    rom_decomp_times, _, _ = read_runtimes(romdecompdirs, "runtime")
    rom_decomp_speedup = [fom_time / time / (domx * domy) for time in rom_decomp_times]
    rom_decomp_speedup_v_decomp = [fom_time_decomp / time for time in rom_decomp_times]
    rom_decomp_errors_plot, rom_decomp_speedup, rom_decomp_speedup_v_decomp = remove_exploded(
        rom_decomp_errors_plot,
        rom_decomp_speedup,
        rom_decomp_speedup_v_decomp
    )

    if plot_hyper:
        hyper_errors, _ = calc_error_norms(
            meshdirs=[meshdir_mono] * (1 + nhyper),
            datalist=fomdata_list + hyperdata_list,
            nvars=nvars[problem],
            dataroot="state_snapshots",
            dtlist=dt[problem],
            samplist=[1] * (nhyper + 1),
            startlist=[startidx_for_rom[problem]] + [0] * nhyper,
            timenorm=True,
            spacenorm=True,
            relative=True,
            mismatch_nan=mismatch_nan,
        )
        hyper_errors_plot = [err[vardata] for err in hyper_errors]
        hyper_times, _, _ = read_runtimes(hyperdirs, "runtime")
        hyper_speedup = [fom_time / time for time in hyper_times]
        hyper_speedup_v_decomp = [fom_time_decomp * domx * domy / time for time in hyper_times]
        hyper_errors_plot, hyper_speedup, hyper_speedup_v_decomp = remove_exploded(
            hyper_errors_plot,
            hyper_speedup,
            hyper_speedup_v_decomp
        )

    if plot_hyper_decomp:
        hyper_decomp_errors, _ = calc_error_norms(
            meshdirs=[meshdir_mono] * (1 + nhyperdecomp),
            datalist=fomdata_list + hyperdecompdata_list,
            nvars=nvars[problem],
            dataroot="state_snapshots",
            dtlist=dt[problem],
            samplist=[1] * (nhyperdecomp + 1),
            startlist=[startidx_for_rom[problem]] + [0] * nhyperdecomp,
            timenorm=True,
            spacenorm=True,
            relative=True,
            mismatch_nan=mismatch_nan,
        )
        hyper_decomp_errors_plot = [err[vardata] for err in hyper_decomp_errors]
        hyperdecomp_times, _, _ = read_runtimes(hyperdecompdirs, "runtime")
        hyperdecomp_speedup = [fom_time / time / (domx * domy) for time in hyperdecomp_times]
        hyperdecomp_speedup_v_decomp = [fom_time_decomp / time for time in hyperdecomp_times]
        hyper_decomp_errors_plot, hyperdecomp_speedup, hyperdecomp_speedup_v_decomp = remove_exploded(
            hyper_decomp_errors_plot,
            hyperdecomp_speedup,
            hyperdecomp_speedup_v_decomp
        )

    plt.sca(ax)
    ax.scatter(rom_speedup_v_decomp, rom_errors_plot, color=plotcolors[0], marker="o")
    if plot_hyper:
        ax.scatter(hyper_speedup_v_decomp, hyper_errors_plot, color=plotcolors[1], marker="o")
    ax.scatter(rom_decomp_speedup_v_decomp, rom_decomp_errors_plot, color=plotcolors[2], marker="o")
    if plot_hyper_decomp:
        ax.scatter(hyperdecomp_speedup_v_decomp, hyper_decomp_errors_plot, color=plotcolors[3], marker="o")

    plt.sca(ax2)
    ax2.scatter(rom_speedup, rom_errors_plot, color=plotcolors[0], marker="o")
    if plot_hyper:
        ax2.scatter(hyper_speedup, hyper_errors_plot, color=plotcolors[1], marker="o")
    ax2.scatter(rom_decomp_speedup, rom_decomp_errors_plot, color=plotcolors[2], marker="o")
    if plot_hyper_decomp:
        ax2.scatter(hyperdecomp_speedup, hyper_decomp_errors_plot, color=plotcolors[3], marker="o")

ax.set_yscale('log')
ax.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xscale('log')

if plot_legend:
    legend_labels = [f"Mono PROM, various M"]
    if plot_hyper:
        legend_labels.append(f"Mono HPROM, fixed M")
    legend_labels.append("Schwarz PROM, various M")
    if plot_hyper_decomp:
        legend_labels.append(f"Schwarz HPROM, fixed M")

    ax.legend(legend_labels, loc="upper left", framealpha=1.0)
    ax2.legend(legend_labels, loc="upper left", framealpha=1.0)


ax.set_xlabel("Speedup vs. Schwarz FOM", labelpad=6.0)
ax.set_ylabel(r"Relative $\ell^2$ error")
ax2.set_xlabel("Speedup vs. Monolithic FOM", labelpad=6.0)
ax.set_xlim(xlim)
ax.set_ylim(error_lims[problem])

if gridlines:
    ax.grid(visible=True, which="major", axis="both", alpha=0.5, color="gray")

plt.tight_layout()
# outfile = os.path.join(outdir, f"pareto_hyper_decomp_overlap_{overlap}_domrate{domrate}_k1_{nmodes_hyper}_k2_{nmodes_hyper_decomp}.png")
outfile = os.path.join(outdir, f"pareto.png")
print(f"Saving image to {outfile}")
plt.savefig(outfile)
print("Finished")
