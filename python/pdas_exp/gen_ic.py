import os

import numpy as np

from pdas.data_utils import load_meshes, decompose_domain_data, write_to_binary
from pdas_exp.defaults import check_params, get_params_combo, nvars
from pdas_exp.utils import mkdir


def gen_ic(
    equations,
    order,
    problem,
    nx,
    ny,
    icFlag,
    icIdx,
    dataroot,
    caseroot,
    phys_params_user={},
    ic_params_user={},
    ndomX=1,
    ndomY=1,
    overlap=None,
):

    numvars = nvars[equations]
    if order == 1:
        order_dir = "firstorder"
    elif order == 3:
        order_dir = "weno3"
    elif order == 5:
        order_dir = "weno5"
    else:
        raise ValueError(f"Invalid order: {order}")

    check_params(phys_params_user, ic_params_user, equations, order, problem, icFlag)
    params_names_list, params_combos = get_params_combo(phys_params_user, ic_params_user, equations, problem, icFlag)

    datadir_base = os.path.join(caseroot, "runs", "fom", f"{nx}x{ny}", order_dir)

    # make directories
    outdir_base = os.path.join(caseroot, "runs", "ic_files")
    mkdir(outdir_base)
    outdir_base = os.path.join(outdir_base, f"{nx}x{ny}")
    mkdir(outdir_base)
    outdir_base = os.path.join(outdir_base, f"{ndomX}x{ndomY}")
    mkdir(outdir_base)

    for param_combo in params_combos:
        # load snapshot
        casedir = ""
        for param_idx, param in enumerate(params_names_list):
            casedir += param + str(param_combo[param_idx]) + "_"
        casedir = casedir[:-1]
        infile = os.path.join(datadir_base, casedir, dataroot + ".bin")
        data = np.fromfile(infile, dtype=np.float64)
        data = np.reshape(data, (numvars, nx, ny, -1), order="F")
        data_snap = data[:, :, :, icIdx]

        if (ndomX == 1) and (ndomY == 1):
            outdir = os.path.join(outdir_base, order_dir)
            mkdir(outdir)
            outdir = os.path.join(outdir, casedir)
            mkdir(outdir)

            data_snap_out = data_snap.flatten(order="F")
            outfile = os.path.join(outdir, f"ic_file_idx{icIdx}.bin")
            print(f"Saving IC file to {outfile}")
            write_to_binary(data_snap_out, outfile)

        else:
            # load mesh
            meshdir = os.path.join(caseroot, "meshes", f"{nx}x{ny}", f"{ndomX}x{ndomY}", f"overlap{overlap}", order_dir, "full")
            _, meshlist_decomp = load_meshes(meshdir, merge_decomp=False)

            # decompose initial conditions
            data_snap_full = np.transpose(data_snap, (1, 2, 0))
            data_decomp = decompose_domain_data(
                data_snap_full,
                meshlist_decomp,
                overlap,
                is_ts=False,
                is_ts_decomp=False,
            )

            outdir = os.path.join(outdir_base, f"overlap{overlap}")
            mkdir(outdir)
            outdir = os.path.join(outdir, order_dir)
            mkdir(outdir)
            outdir = os.path.join(outdir, casedir)
            mkdir(outdir)

            # write to binary
            print(f"Saving IC file to {outdir}")
            for j in range(ndomY):
                for i in range(ndomX):
                    dom_idx = i + 2 * j
                    data_out = np.transpose(data_decomp[i][j][0], (2, 0, 1))
                    data_out = data_out.flatten(order="F")
                    outfile = os.path.join(outdir, f"ic_file_idx{icIdx}_dom{dom_idx}.bin")
                    write_to_binary(data_out, outfile)
