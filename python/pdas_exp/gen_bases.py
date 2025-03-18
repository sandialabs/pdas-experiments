
import os
import argparse

import yaml

from pschwarz.prom_utils import gen_pod_bases
from pdas_exp.defaults import check_params, get_params_combo, nvars
from pdas_exp.utils import check_meshdir, mkdir


def gen_bases(
    equations,
    order,
    problem,
    meshroot,
    nx,
    ny,
    icFlag,
    rundir_base,
    dataroot,
    center_method,
    norm_method,
    nmodes,
    outdir,
    phys_params_user={},
    ic_params_user={},
    idx_start=0,
    idx_end=None,
    idx_skip=1,
    decomp=True,
    ndomX=-1,
    ndomY=-1,
    overlap=-1,
    recalc=True,
):

    check_params(phys_params_user, ic_params_user, equations, order, problem, icFlag)

    if order == 1:
        order_dir = "firstorder"
    elif order == 3:
        order_dir = "weno3"
    elif order == 5:
        order_dir = "weno5"
    else:
        raise ValueError(f"Invalid order: {order}")

    # generate output directories
    mkdir(outdir)
    outdir = os.path.join(outdir, f"{nx}x{ny}")
    mkdir(outdir)
    if decomp:
        assert ndomX >= 1
        assert ndomY >= 1
        assert overlap >= 0
        ndomains = ndomX * ndomY
        assert ndomains >= 2
        outdir = os.path.join(outdir, f"{ndomX}x{ndomY}")
        mkdir(outdir)
        outdir = os.path.join(outdir, f"overlap{overlap}")
        mkdir(outdir)
    else:
        outdir = os.path.join(outdir, "1x1")
        mkdir(outdir)

    # before final directory, if it already exists don't redo
    outdir = os.path.join(outdir, order_dir)
    if os.path.isdir(outdir) and (not recalc):
        print(f"Not recalculating basis at {outdir}")
        return
    else:
        mkdir(outdir)

    # check mesh(es)
    meshdir = os.path.join(meshroot, f"{nx}x{ny}", "1x1", order_dir, "full")
    assert os.path.isdir(meshdir)
    check_meshdir(meshdir)
    if decomp:
        meshdir_decomp = os.path.join(meshroot, f"{nx}x{ny}", f"{ndomX}x{ndomY}", f"overlap{overlap}", order_dir, "full")
        assert os.path.isfile(os.path.join(meshdir_decomp, "info_domain.dat"))
        for dom_idx in range(ndomains):
            check_meshdir(os.path.join(meshdir_decomp, f"domain_{dom_idx}"))
    else:
        meshdir_decomp = None

    # gather data file paths
    params_names_list, params_combo = get_params_combo(phys_params_user, ic_params_user, equations, problem, icFlag)
    nruns = len(params_combo)
    datadirs = []
    for run_idx, run_list in enumerate(params_combo):
        rundir = os.path.join(rundir_base, "fom", f"{nx}x{ny}", order_dir)
        dirname = ""
        for param_idx, param in enumerate(params_names_list):
            dirname += param + str(run_list[param_idx]) + "_"
        dirname = dirname[:-1]
        datadirs.append(os.path.join(rundir, dirname))

    gen_pod_bases(
        outdir,
        meshdir=[meshdir]*nruns,
        datadir=datadirs,
        nvars=nvars[equations],
        dataroot=dataroot,
        concat=True,
        pod_decomp=decomp,
        meshdir_decomp=meshdir_decomp,
        idx_start=idx_start,
        idx_end=idx_end,
        idx_skip=idx_skip,
        center_method=center_method,
        norm_method=norm_method,
        nmodes=nmodes,
    )

    print("Finished")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    args = parser.parse_args()
    f = open(args.settings_file, "r")
    inputs = yaml.safe_load(f)

    # handle parameters
    phys_params_user = inputs["phys_params_user"]
    if phys_params_user is None:
        phys_params_user = {}
    elif isinstance(phys_params_user, dict):
        for key, value in phys_params_user.items():
            if not isinstance(value, list):
                phys_params_user[key] = [value]
    else:
        raise ValueError("phys_params_user must be a nested input")

    ic_params_user = inputs["ic_params_user"]
    if ic_params_user is None:
        ic_params_user = {}
    elif isinstance(ic_params_user, dict):
        for key, value in ic_params_user.items():
            if not isinstance(value, list):
                ic_params_user[key] = [value]
    else:
        raise ValueError("ic_params_user must be a nested input")

    # handle sampling indices
    try:
        idx_start = inputs["idx_start"]
    except KeyError:
        idx_start = 1
    try:
        idx_end = inputs["idx_end"]
    except KeyError:
        idx_end = None
    try:
        idx_skip = inputs["idx_skip"]
    except KeyError:
        idx_skip = 1

    if inputs["decomp"]:
        ndomX = inputs["ndomX"]
        ndomY = inputs["ndomY"]
        overlap = inputs["overlap"]
    else:
        ndomX = None
        ndomY = None
        overlap = None

    gen_bases(
        inputs["equations"],
        inputs["order"],
        inputs["problem"],
        inputs["meshroot"],
        inputs["nx"],
        inputs["ny"],
        inputs["icFlag"],
        inputs["rundir_base"],
        inputs["dataroot"],
        inputs["center_method"],
        inputs["norm_method"],
        inputs["nmodes"],
        inputs["outdir"],
        phys_params_user=phys_params_user,
        ic_params_user=ic_params_user,
        idx_start=idx_start,
        idx_end=idx_end,
        idx_skip=idx_skip,
        decomp=inputs["decomp"],
        ndomX=ndomX,
        ndomY=ndomY,
        overlap=overlap,
    )