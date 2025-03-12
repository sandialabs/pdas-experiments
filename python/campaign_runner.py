import os
import argparse
import subprocess

import yaml
import numpy as np

from pschwarz.samp_utils import gen_sample_mesh
from pdas_exp.utils import catchinput, mkdir
from pdas_exp.gen_runs import gen_runs
from pdas_exp.gen_bases import gen_bases
from pdas_exp.gen_ic import gen_ic


def main(inputs):

    # read campaign configuration parameters
    hotrun = catchinput(inputs, "hotrun", default_val=True)
    recalc = catchinput(inputs, "recalc", default_val=False)
    runexec_serial = catchinput(inputs, "runexec_serial", default_type=str)
    assert os.path.isfile(runexec_serial)
    runexec_omp = catchinput(inputs, "runexec_omp", default_type=str)
    assert os.path.isfile(runexec_omp)
    pdaroot   = catchinput(inputs, "pdaroot", default_type=str)
    assert os.path.isdir(pdaroot)
    pdasroot  = catchinput(inputs, "pdasroot", default_type=str)
    assert os.path.isdir(pdasroot)
    caseroot   = catchinput(inputs, "caseroot", default_type=str)
    assert os.path.isdir(caseroot)

    # read global problem inputs
    inputs_case = inputs["case"]
    equations     = catchinput(inputs_case, "equations",     default_type=str)
    problem       = catchinput(inputs_case, "problem",          default_type=str)
    ic_flag       = catchinput(inputs_case, "ic_flag",       default_type=int)
    flux_order    = catchinput(inputs_case, "flux_order",    default_type=int)
    time_scheme   = catchinput(inputs_case, "time_scheme",   default_type=str)
    tf            = catchinput(inputs_case, "tf",            default_type=float)
    dt            = catchinput(inputs_case, "dt",            default_type=float)

    # read mesh inputs
    inputs_mesh = inputs["mesh"]
    nx      = catchinput(inputs_mesh, "nx",      default_type=int)
    ny      = catchinput(inputs_mesh, "ny",      default_type=int)
    xbounds = catchinput(inputs_mesh, "xbounds", default_type=list)
    ybounds = catchinput(inputs_mesh, "ybounds", default_type=list)

    # parameterizations
    ic_parameterized = "ic_param" in inputs
    ic_params = {}
    if ic_parameterized:
        inputs_ic_param = inputs["ic_param"]
        ic_param_name = catchinput(inputs_ic_param, "name", default_type=str)
        ic_param_vals = catchinput(inputs_ic_param, "vals", default_type=list)
        ic_params = {ic_param_name: ic_param_vals}

    phys_parameterized = "phys_param" in inputs
    phys_params = {}
    if phys_parameterized:
        inputs_phys_param = inputs["phys_param"]
        phys_param_name = catchinput(inputs_phys_param, "name", default_type=str)
        phys_param_vals = catchinput(inputs_phys_param, "vals", default_type=list)
        phys_params = {phys_param_name: phys_param_vals}

    assert not (ic_parameterized and phys_parameterized)

    # run flags and inputs
    run_fom = "fom" in inputs
    if run_fom:
        inputs_fom = inputs["fom"]
    run_prom = "prom" in inputs
    if run_prom:
        inputs_prom = inputs["prom"]
    run_hprom = "hprom" in inputs
    if run_hprom:
        inputs_hprom = inputs["hprom"]
    run_fom_decomp = "fom_decomp" in inputs
    if run_fom_decomp:
        inputs_fom_decomp = inputs["fom_decomp"]
    run_prom_decomp  = "prom_decomp" in inputs
    if run_prom_decomp:
        inputs_prom_decomp = inputs["prom_decomp"]
    run_hprom_decomp = "hprom_decomp" in inputs
    if run_hprom_decomp:
        inputs_hprom_decomp = inputs["hprom_decomp"]

    # meshing
    mesh_driver = os.path.join(pdaroot, "meshing_scripts", "create_full_mesh.py")
    assert os.path.isfile(mesh_driver)
    decomp_driver = os.path.join(pdasroot, "meshing_scripts", "create_decomp_meshes.py")
    assert os.path.isfile(decomp_driver)
    sample_driver = os.path.join(pdaroot, "meshing_scripts", "create_sample_mesh.py")
    assert os.path.isfile(sample_driver)
    meshroot = os.path.join(caseroot, "meshes")
    mkdir(meshroot)
    if flux_order == 1:
        stencildir = "firstorder"
        stencil = 3
    elif flux_order == 3:
        stencildir = "weno3"
        stencil = 5
    elif flux_order == 5:
        stencildir = "weno5"
        stencil = 7
    else:
        raise ValueError(f"Invalid flux_order: {flux_order}")

    # monolithic mesh
    if run_fom or run_prom or run_hprom:
        meshpath = os.path.join(meshroot, f"{nx}x{ny}", "1x1", stencildir, "full")
        if not os.path.isdir(meshpath):
            os.makedirs(meshpath, exist_ok=True)
            print(f"Generating mesh at {meshpath}")
            mesh_cmd = f"python3 {mesh_driver} "
            mesh_cmd += f"--numCells {nx} {ny} "
            mesh_cmd += f"--outDir {meshpath} "
            mesh_cmd += f"-s {stencil} "
            mesh_cmd += f"--bounds {xbounds[0]} {xbounds[1]} {ybounds[0]} {ybounds[1]}"
            subprocess.run([mesh_cmd], shell=True)

    # decomposed mesh
    # determine all necessary overlap values
    overlaps_mesh = []
    if run_fom_decomp:
        overlaps_fom_decomp = catchinput(inputs_fom_decomp, "overlaps", default_type=list)
        overlaps_mesh = np.union1d(overlaps_mesh, overlaps_fom_decomp)
    if run_prom_decomp:
        overlaps_prom_decomp = catchinput(inputs_prom_decomp, "overlaps", default_type=list)
        overlaps_mesh = np.union1d(overlaps_mesh, overlaps_prom_decomp)
    if run_hprom_decomp:
        overlaps_hprom_decomp = catchinput(inputs_hprom_decomp, "overlaps", default_type=list)
        overlaps_mesh = np.union1d(overlaps_mesh, overlaps_hprom_decomp)
    overlaps_mesh = [round(overlap) for overlap in overlaps_mesh]
    if run_fom_decomp or run_prom_decomp or run_hprom_decomp:
        # Schwarz parameters
        inputs_schwarz = inputs["schwarz"]
        domx          = catchinput(inputs_schwarz, "domx", default_type=int)
        domy          = catchinput(inputs_schwarz, "domx", default_type=int)
        isadditive    = catchinput(inputs_schwarz, "isadditive", default_type=bool)
        conv_step_max = catchinput(inputs_schwarz, "conv_step_max", default_type=int)

        numprocs_omp = domx * domy

        for overlap in overlaps_mesh:
            meshpath = os.path.join(meshroot, f"{nx}x{ny}", f"{domx}x{domy}", f"overlap{overlap}", stencildir, "full")
            if not os.path.isdir(meshpath):
                os.makedirs(meshpath, exist_ok=True)
                print(f"Generating mesh at {meshpath}")
                mesh_cmd = f"python3 {decomp_driver} "
                mesh_cmd += f"--meshscript {mesh_driver} "
                mesh_cmd += f"--outdir {meshpath} "
                mesh_cmd += f"--numcells {nx} {ny} "
                mesh_cmd += f"--bounds {xbounds[0]} {xbounds[1]} {ybounds[0]} {ybounds[1]} "
                mesh_cmd += f"--stencilsize {stencil} "
                mesh_cmd += f"--numdoms {domx} {domy} "
                mesh_cmd += f"--overlap {overlap}"
                subprocess.run([mesh_cmd], shell=True)

    runroot = os.path.join(caseroot, "runs")
    mkdir(runroot)

    # run FOMs
    # monolithic only iterate over parameterization
    if run_fom:
        sampfreq = catchinput(inputs_fom, "sampfreq", default_type=int)
        gen_runs(
            equations,
            flux_order,
            problem,
            time_scheme,
            tf,
            dt,
            sampfreq,
            meshroot,
            nx,
            ny,
            "fom",
            ic_flag,
            runexec_serial,
            "info",
            "file",
            runroot,
            phys_params_user=phys_params,
            ic_params_user=ic_params,
            run=hotrun,
            recalc=recalc,
        )

    # decomp iterates over overlap
    if run_fom_decomp:
        sampfreq = catchinput(inputs_fom_decomp, "sampfreq", default_type=int)
        overlaps = catchinput(inputs_fom_decomp, "overlaps", default_type=list)
        for overlap in overlaps:
            gen_runs(
                equations,
                flux_order,
                problem,
                time_scheme,
                tf,
                dt,
                sampfreq,
                meshroot,
                nx,
                ny,
                "fom_decomp",
                ic_flag,
                runexec_omp,
                "info",
                "file",
                runroot,
                phys_params_user=phys_params,
                ic_params_user=ic_params,
                solve_algo="FOM",
                ndomX=domx,
                ndomY=domy,
                overlap=overlap,
                isadditive=isadditive,
                conv_step_max=conv_step_max,
                numprocs=numprocs_omp,
                run=hotrun,
                recalc=recalc,
            )

    # generate ROM ICs, if requested
    if "rom_from_ic" in inputs:
        ic_index = catchinput(inputs["rom_from_ic"], "index", default_type=int)

        # account for final time change
        tf_rom = tf - ic_index * dt

        # monolithic
        gen_ic(
            equations,
            flux_order,
            problem,
            nx,
            ny,
            ic_flag,
            ic_index,
            "state_snapshots",
            caseroot,
            phys_params_user=phys_params,
            ic_params_user=ic_params,
        )
        # decomposed
        for overlap in overlaps_mesh:
            gen_ic(
                equations,
                flux_order,
                problem,
                nx,
                ny,
                ic_flag,
                ic_index,
                "state_snapshots",
                caseroot,
                phys_params_user=phys_params,
                ic_params_user=ic_params,
                ndomX=domx,
                ndomY=domy,
                overlap=overlap,
            )
    else:
        ic_index = None
        tf_rom = tf

    # generate trial bases
    input_basis = inputs["basis"]
    basis_name    = catchinput(input_basis, "name", default_type=str)
    basis_vals    = catchinput(input_basis, "param_vals", default_type=list)
    nmodes_max    = catchinput(input_basis, "nmodes_max", default_type=int)
    idx_start     = catchinput(input_basis, "idx_start", default_type=int)
    idx_end       = catchinput(input_basis, "idx_end", default_type=int)
    idx_skip      = catchinput(input_basis, "idx_skip", default_type=int)
    center_method = catchinput(input_basis, "center_method", default_type=str)
    norm_method   = catchinput(input_basis, "norm_method", default_type=str)

    basisroot = os.path.join(caseroot, "pod_bases")
    mkdir(basisroot)
    basisdir_trial = os.path.join(basisroot, basis_name)
    mkdir(basisdir_trial)

    ic_params_basis = {}
    if len(ic_params) != 0:
        ic_params_basis = {ic_param_name: basis_vals}
    phys_params_basis = {}
    if len(phys_params) != 0:
        phys_params_basis = {phys_param_name: basis_vals}

    # monolithic bases
    if run_prom or run_hprom:
        gen_bases(
            equations,
            flux_order,
            problem,
            meshroot,
            nx,
            ny,
            ic_flag,
            runroot,
            "state_snapshots",
            center_method,
            norm_method,
            nmodes_max,
            basisdir_trial,
            phys_params_user=phys_params_basis,
            ic_params_user=ic_params_basis,
            idx_start=idx_start,
            idx_end=idx_end,
            idx_skip=idx_skip,
            decomp=False,
            recalc=recalc,
        )

    # decomposed bases
    if run_prom_decomp or run_hprom_decomp:
        for overlap in overlaps_mesh:
            gen_bases(
                equations,
                flux_order,
                problem,
                meshroot,
                nx,
                ny,
                ic_flag,
                runroot,
                "state_snapshots",
                center_method,
                norm_method,
                nmodes_max,
                basisdir_trial,
                phys_params_user=phys_params_basis,
                ic_params_user=ic_params_basis,
                idx_start=idx_start,
                idx_end=idx_end,
                idx_skip=idx_skip,
                decomp=True,
                ndomX=domx,
                ndomY=domy,
                overlap=overlap,
                recalc=recalc,
            )

    # PROM runs
    # monolithic iterates over nmodes
    if run_prom:
        rom_algo     = catchinput(inputs_prom, "rom_algo", default_type=str)
        sampfreq     = catchinput(inputs_prom, "sampfreq", default_type=int)
        nmodes_trial = catchinput(inputs_prom, "nmodes_trial", default_type=list)
        for nmodes in nmodes_trial:
            gen_runs(
                equations,
                flux_order,
                problem,
                time_scheme,
                tf_rom,
                dt,
                sampfreq,
                meshroot,
                nx,
                ny,
                "rom",
                ic_flag,
                runexec_serial,
                "info",
                "file",
                runroot,
                phys_params_user=phys_params,
                ic_params_user=ic_params,
                solve_algo=rom_algo,
                nmodes=nmodes,
                basis_dir=basisdir_trial,
                basis_file="basis",
                shift_file="center",
                ic_index=ic_index,
                run=hotrun,
                recalc=recalc,
            )

    # decomposed iterates over modes and overlap
    if run_prom_decomp:
        rom_algo     = catchinput(inputs_prom_decomp, "rom_algo", default_type=str)
        sampfreq     = catchinput(inputs_prom_decomp, "sampfreq", default_type=int)
        nmodes_trial = catchinput(inputs_prom_decomp, "nmodes_trial", default_type=list)
        overlaps     = catchinput(inputs_prom_decomp, "overlaps", default_type=list)
        for nmodes in nmodes_trial:
            for overlap in overlaps:
                gen_runs(
                    equations,
                    flux_order,
                    problem,
                    time_scheme,
                    tf_rom,
                    dt,
                    sampfreq,
                    meshroot,
                    nx,
                    ny,
                    "rom_decomp",
                    ic_flag,
                    runexec_omp,
                    "info",
                    "file",
                    runroot,
                    phys_params_user=phys_params,
                    ic_params_user=ic_params,
                    solve_algo=rom_algo,
                    nmodes=nmodes,
                    basis_dir=basisdir_trial,
                    basis_file="basis",
                    shift_file="center",
                    ndomX=domx,
                    ndomY=domy,
                    overlap=overlap,
                    isadditive=isadditive,
                    conv_step_max=conv_step_max,
                    numprocs=numprocs_omp,
                    ic_index=ic_index,
                    run=hotrun,
                    recalc=recalc,
                )

    # monolithic HPROM
    if run_hprom:
        # solver inputs
        rom_algo = catchinput(inputs_hprom, "rom_algo", default_type=str)
        nmodes_trial = catchinput(inputs_hprom, "nmodes_trial", default_type=list)
        sampfreq = catchinput(inputs_hprom, "sampfreq", default_type=int)
        gpod_weigher = catchinput(inputs_hprom, "gpod_weigher", default_type=str)
        if gpod_weigher != "identity":
            basis_name_gpod = catchinput(inputs_hprom, "basis_name_gpod", default_type=str)
            basisdir_gpod = os.path.join(basisroot, basis_name_gpod)
            nmodes_gpod = catchinput(inputs_hprom, "nmodes_gpod", default_type=int)
        else:
            basisdir_gpod = None
            nmodes_gpod = None

        # sample mesh inputs
        sampalgo = catchinput(inputs_hprom, "sampalgo", default_type=str)
        qdeim = catchinput(inputs_hprom, "qdeim", default_type=bool)
        seed_phys = catchinput(inputs_hprom, "seed_phys", default_type=bool)
        if seed_phys:
            seed_phys_rate = catchinput(inputs_hprom, "seed_phys_rate", default_type=int)
        else:
            seed_phys_rate = 0
        samp_phys = catchinput(inputs_hprom, "samp_phys", default_type=bool)
        samp_percs = catchinput(inputs_hprom, "samp_percs", default_type=list)

        meshdir = os.path.join(
            meshroot,
            f"{nx}x{ny}",
            "1x1",
            stencildir,
        )
        meshdir_full = os.path.join(meshdir, "full")
        meshroot_samp = os.path.join(meshdir, sampalgo)
        mkdir(meshroot_samp)

        if qdeim or (sampalgo in ["eigenvec", "gnat"]):
            raise ValueError("Eigenvec not implemented yet")
            # add too meshroot_samp
        else:
            basis_dir_greedy = None
            nmodes_greedy = None
            meshroot_samp += "/"

        for samp_perc in samp_percs:
            # generate sample mesh
            outdir = meshroot_samp + f"samp_{samp_perc}"
            if qdeim:
                outdir += "_qdeim"
            if seed_phys:
                outdir += f"_phys{seed_phys_rate}"

            gen_sample_mesh(
                sampalgo,
                meshdir_full,
                samp_perc,
                outdir,
                basis_dir=basis_dir_greedy,
                nmodes=nmodes_greedy,
                seed_qdeim=qdeim,
                seed_phys_bounds=seed_phys,
                seed_phys_rate=seed_phys_rate,
                seed_dom_bounds=False,
                seed_dom_rate=0,
                samp_phys_bounds=samp_phys,
                samp_dom_bounds=False,
                recalc=recalc,
            )
            print(f"Generating mesh at {meshpath}")
            mesh_cmd = f"python3 {sample_driver} "
            mesh_cmd += f"--fullMeshDir {meshdir_full} "
            mesh_cmd += f"--sampleMeshIndices {outdir}/sample_mesh_gids.dat "
            mesh_cmd += f"--outDir {outdir}"
            subprocess.run([mesh_cmd], shell=True)

            # run HPROM, already iterating over samp_perc
            # can iterate over trial modes
            for nmodes in nmodes_trial:
                gen_runs(
                    equations,
                    flux_order,
                    problem,
                    time_scheme,
                    tf_rom,
                    dt,
                    sampfreq,
                    meshroot,
                    nx,
                    ny,
                    "hyper",
                    ic_flag,
                    runexec_serial,
                    "info",
                    "file",
                    runroot,
                    phys_params_user=phys_params,
                    ic_params_user=ic_params,
                    solve_algo=rom_algo,
                    nmodes=nmodes,
                    basis_dir=basisdir_trial,
                    basis_file="basis",
                    shift_file="center",
                    ic_index=ic_index,
                    sampalgo=sampalgo,
                    nmodes_greedy=nmodes_greedy,
                    gpod_weigher=gpod_weigher,
                    basis_dir_gpod=basisdir_gpod,
                    nmodes_gpod=nmodes_gpod,
                    seed_qdeim=qdeim,
                    sampperc=samp_perc,
                    physrate=seed_phys_rate,
                    run=hotrun,
                    recalc=recalc,
                )

    # monolithic HPROM
    if run_hprom_decomp:
        # solver inputs
        rom_algo = catchinput(inputs_hprom_decomp, "rom_algo", default_type=str)
        nmodes_trial = catchinput(inputs_hprom_decomp, "nmodes_trial", default_type=list)
        sampfreq = catchinput(inputs_hprom_decomp, "sampfreq", default_type=int)
        gpod_weigher = catchinput(inputs_hprom_decomp, "gpod_weigher", default_type=str)
        if gpod_weigher != "identity":
            basis_name_gpod = catchinput(inputs_hprom_decomp, "basis_name_gpod", default_type=str)
            basisdir_gpod = os.path.join(basisroot, basis_name_gpod)
            nmodes_gpod = catchinput(inputs_hprom_decomp, "nmodes_gpod", default_type=int)
        else:
            basisdir_gpod = None
            nmodes_gpod = None

        # sample mesh inputs
        sampalgo = catchinput(inputs_hprom_decomp, "sampalgo", default_type=str)
        qdeim = catchinput(inputs_hprom_decomp, "qdeim", default_type=bool)
        seed_phys = catchinput(inputs_hprom_decomp, "seed_phys", default_type=bool)
        if seed_phys:
            seed_phys_rate = catchinput(inputs_hprom_decomp, "seed_phys_rate", default_type=int)
        else:
            seed_phys_rate = 0
        samp_phys = catchinput(inputs_hprom_decomp, "samp_phys", default_type=bool)
        seed_dom = catchinput(inputs_hprom_decomp, "seed_dom", default_type=bool)
        if seed_dom:
            seed_dom_rates = catchinput(inputs_hprom_decomp, "seed_dom_rates", default_type=list)
        else:
            seed_dom_rate = [0]
        samp_dom = catchinput(inputs_hprom_decomp, "samp_dom", default_type=bool)
        samp_percs = catchinput(inputs_hprom_decomp, "samp_percs", default_type=list)

        # can iterate over overlap, samp_perc, and domrate
        for overlap in overlaps_hprom_decomp:

            meshdir = os.path.join(
                meshroot,
                f"{nx}x{ny}",
                f"{domx}x{domy}",
                f"overlap{overlap}",
                stencildir,
            )
            meshdir_full = os.path.join(meshdir, "full")
            meshroot_samp = os.path.join(meshdir, sampalgo)
            mkdir(meshroot_samp)

            if qdeim or (sampalgo in ["eigenvec", "gnat"]):
                raise ValueError("Eigenvec not implemented yet")
                # add too meshroot_samp
            else:
                basis_dir_greedy = None
                nmodes_greedy = None
                meshroot_samp += "/"

            for samp_perc in samp_percs:
                for seed_dom_rate in seed_dom_rates:
                    # generate sample mesh
                    outdir = meshroot_samp + f"samp_{samp_perc}"
                    if qdeim:
                        outdir += "_qdeim"
                    if seed_phys:
                        outdir += f"_phys{seed_phys_rate}"
                    if seed_dom:
                        outdir += f"_dom{seed_dom_rate}"
                    mkdir(outdir)

                    gen_sample_mesh(
                        sampalgo,
                        meshdir_full,
                        samp_perc,
                        outdir,
                        basis_dir=basis_dir_greedy,
                        nmodes=nmodes_greedy,
                        seed_qdeim=qdeim,
                        seed_phys_bounds=seed_phys,
                        seed_phys_rate=seed_phys_rate,
                        seed_dom_bounds=seed_dom,
                        seed_dom_rate=seed_dom_rate,
                        samp_phys_bounds=samp_phys,
                        samp_dom_bounds=samp_dom,
                        recalc=recalc,
                    )

                    # run decomposed HPROM
                    # further iterate over trial modes
                    for nmodes in nmodes_trial:
                        gen_runs(
                            equations,
                            flux_order,
                            problem,
                            time_scheme,
                            tf_rom,
                            dt,
                            sampfreq,
                            meshroot,
                            nx,
                            ny,
                            "hyper_decomp",
                            ic_flag,
                            runexec_omp,
                            "info",
                            "file",
                            runroot,
                            phys_params_user=phys_params,
                            ic_params_user=ic_params,
                            solve_algo=rom_algo,
                            nmodes=nmodes,
                            basis_dir=basisdir_trial,
                            basis_file="basis",
                            shift_file="center",
                            ndomX=domx,
                            ndomY=domy,
                            overlap=overlap,
                            isadditive=isadditive,
                            conv_step_max=conv_step_max,
                            numprocs=numprocs_omp,
                            ic_index=ic_index,
                            sampalgo=sampalgo,
                            nmodes_greedy=nmodes_greedy,
                            gpod_weigher=gpod_weigher,
                            basis_dir_gpod=basisdir_gpod,
                            nmodes_gpod=nmodes_gpod,
                            seed_qdeim=qdeim,
                            sampperc=samp_perc,
                            physrate=seed_phys_rate,
                            domrate=seed_dom_rate,
                            run=hotrun,
                            recalc=recalc,
                        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    args = parser.parse_args()

    with open(args.settings_file, "r") as f:
        inputs = yaml.safe_load(f)

    main(inputs)