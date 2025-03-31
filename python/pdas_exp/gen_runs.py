
import os
import subprocess
import argparse

import yaml

from .defaults import check_params, get_params_combo
from .utils import check_meshdir, mkdir, catchlist

ALGOS = ["FOM", "LSPG", "LSPGHyper"]

def gen_runs(
    equations,
    order,
    problem,
    scheme,
    tf,
    dt,
    sampfreq,
    meshroot,
    nx,
    ny,
    runtype,
    icFlag,
    runner,
    loglevel,
    logtarget,
    outdir_base,
    phys_params_user={},
    ic_params_user={},
    solve_algo=None,
    nmodes=None,
    basis_dir=None,
    basis_file=None,
    shift_file=None,
    ndomX=None,
    ndomY=None,
    overlap=None,
    isadditive=False,
    conv_step_max=10,
    numprocs=1,
    ic_index=None,
    sampalgo=None,
    nmodes_greedy=None,
    gpod_weigher=None,
    basis_dir_gpod=None,
    nmodes_gpod=None,
    seed_qdeim=False,
    sampperc=0.0,
    physrate=0,
    domrate=0,
    run=False,
    recalc=True,
):

    # ----- START CHECKS -----

    assert loglevel in ["debug", "info", "none"]
    assert logtarget in ["file", "console", "both"]
    assert scheme in ["BDF1", "BDF2"]
    assert tf > 0.0
    assert sampfreq >= 1
    assert runtype in [
        "fom", "rom", "hyper",
        "fom_decomp", "rom_decomp", "hyper_decomp",
    ]
    assert os.path.isdir(outdir_base)
    assert os.path.isdir(meshroot)
    if ic_index is not None:
        assert ic_index >= 0

    check_params(phys_params_user, ic_params_user, equations, order, problem, icFlag)

    if "hyper" in runtype:
        assert gpod_weigher in ["identity", "gappy_pod"]
        assert sampalgo in ["random", "eigenvec"]
        assert sampperc > 0.0

    if order == 1:
        order_dir = "firstorder"
    elif order == 3:
        order_dir = "weno3"
    elif order == 5:
        order_dir = "weno5"
    else:
        raise ValueError(f"Invalid order: {order}")

    if "decomp" not in runtype:
        ndomX = 1
        ndomY = 1
        ndomains = 1
        assert dt > 0.0
    else:
        assert ndomX >= 1
        assert ndomY >= 1
        ndomains = ndomX * ndomY
        assert ndomains > 1, "No point to decomp if 1x1"
        assert overlap is not None
        assert overlap >= 0
        assert numprocs > 0
        if not isadditive:
            assert numprocs == 1
        os.environ["OMP_NUM_THREADS"] = str(numprocs)
        dt = catchlist(dt, float, ndomains)

    # handle nmodes, algorithm
    if runtype != "fom":
        assert solve_algo is not None

        if runtype != "fom_decomp":
            assert nmodes is not None
            assert basis_dir is not None
            assert basis_file is not None
            assert shift_file is not None

        if "decomp" in runtype:
            solve_algo = catchlist(solve_algo, str, ndomains)
            assert all([algo in ALGOS for algo in solve_algo])
            if runtype == "fom_decomp":
                assert all([algo == "FOM" for algo in solve_algo])
            else:
                nmodes = catchlist(nmodes, int, ndomains)
                basis_dir = os.path.join(basis_dir, f"{nx}x{ny}", f"{ndomX}x{ndomY}")
                basis_dir = os.path.join(basis_dir, f"overlap{overlap}", order_dir)
                basis_root = os.path.join(basis_dir, basis_file)
                shift_root = os.path.join(basis_dir, shift_file)
                assert all([
                    os.path.isfile(f"{basis_root}_{dom_idx}.bin") \
                    for dom_idx in range(ndomains)
                ]), f"Basis file not found at {basis_dir}"
                assert all([
                    os.path.isfile(f"{shift_root}_{dom_idx}.bin") \
                    for dom_idx in range(ndomains)
                ]), f"Affine shift file not found at {basis_dir}"
        else:
            assert isinstance(nmodes, int)
            assert solve_algo in ALGOS
            basis_dir = os.path.join(basis_dir, f"{nx}x{ny}", "1x1", order_dir)
            basis_root = os.path.join(basis_dir, basis_file)
            shift_root = os.path.join(basis_dir, shift_file)
            assert os.path.isfile(basis_root + ".bin"), f"No basis at {basis_root}.bin"
            assert os.path.isfile(shift_root + ".bin"), f"No shift at {shift_root}.bin"

    # check that mesh exists
    meshdir = os.path.join(meshroot, f"{nx}x{ny}", f"{ndomX}x{ndomY}")
    assert os.path.isdir(meshdir)

    if "decomp" not in runtype:
        meshdir_full = os.path.join(meshdir, order_dir, "full")
        check_meshdir(meshdir_full)
    else:
        meshdir_full = os.path.join(meshdir, f"overlap{overlap}", order_dir, "full")
        assert os.path.isfile(os.path.join(meshdir_full, "info_domain.dat"))
        for dom_idx in range(ndomains):
            check_meshdir(os.path.join(meshdir_full, f"domain_{dom_idx}"))

    if "hyper" in runtype:

        if runtype == "hyper":
            finaldir = ""
            if (sampalgo != "random") or seed_qdeim:
                assert isinstance(nmodes_greedy, int)
                finaldir += f"modes_{nmodes_greedy}_"
            finaldir += f"samp_{sampperc}"
            if seed_qdeim:
                finaldir += "_qdeim"
            meshdir_hyper = os.path.join(meshdir, order_dir, sampalgo, finaldir)
            if gpod_weigher != "identity":
                basis_dir_gpod = os.path.join(basis_dir_gpod, f"{nx}x{ny}", "1x1", order_dir)
                basis_root_gpod = os.path.join(basis_dir_gpod, basis_file)

        elif runtype == "hyper_decomp":
            finaldir = ""
            if (sampalgo != "random") or seed_qdeim:
                nmodes_greedy = catchlist(nmodes_greedy, int, ndomains)
                assert all([modeval == nmodes_greedy[0] for modeval in nmodes_greedy]), "Can't handle different modes in each hyp-red domain yet"
                finaldir += f"modes_{nmodes_greedy[0]}_"
            finaldir += f"samp_{sampperc}"
            if seed_qdeim:
                finaldir += "_qdeim"
            meshdir_hyper = os.path.join(meshdir, f"overlap{overlap}", order_dir, sampalgo, finaldir)
            if gpod_weigher != "identity":
                basis_dir_gpod = os.path.join(basis_dir_gpod, f"{nx}x{ny}", f"{ndomX}x{ndomY}")
                basis_dir_gpod = os.path.join(basis_dir_gpod, f"overlap{overlap}", order_dir)
                basis_root_gpod = os.path.join(basis_dir_gpod, basis_file)

        if physrate > 0:
            meshdir_hyper += f"_phys{physrate}"
        if domrate > 0:
            meshdir_hyper += f"_dom{domrate}"
        assert os.path.isdir(meshdir_hyper), meshdir_hyper

        if runtype == "hyper_decomp":
            sampfiles = [os.path.join(meshdir_hyper, f"domain_{dom_idx}", "sample_mesh_gids.dat") for dom_idx in range(ndomains)]

    # ----- END CHECKS -----

    # ----- START RUN GENERATION -----

    # generate permutations of parameter lists
    params_names_list, params_combo = get_params_combo(phys_params_user, ic_params_user, equations, problem, icFlag)

    for run_idx, run_list in enumerate(params_combo):

        rundir = os.path.join(outdir_base, runtype)
        mkdir(rundir)

        # mesh directory
        rundir = os.path.join(rundir, f"{nx}x{ny}")
        mkdir(rundir)

        # decomp directory
        if "decomp" in runtype:
            rundir = os.path.join(rundir, f"{ndomX}x{ndomY}")
            mkdir(rundir)
            rundir = os.path.join(rundir, f"overlap{overlap}")
            mkdir(rundir)

        rundir = os.path.join(rundir, order_dir)
        mkdir(rundir)

        # parameter directory
        paramdir = ""
        for param_idx, param in enumerate(params_names_list):
            paramdir += f"{param}{run_list[param_idx]}_"
        paramdir = paramdir[:-1]
        rundir = os.path.join(rundir, paramdir)

        # find initial conditions file, if requested
        if ic_index is not None:
            ic_dir = os.path.join(outdir_base, "ic_files", f"{nx}x{ny}", f"{ndomX}x{ndomY}")
            if "decomp" in runtype:
                ic_dir  = os.path.join(ic_dir, f"overlap{overlap}", order_dir, paramdir)
                ic_file = f"ic_file_idx{ic_index}_dom"
            else:
                ic_dir = os.path.join(ic_dir, order_dir, paramdir)
                ic_file = f"ic_file_idx{ic_index}.bin"
            ic_file = os.path.join(ic_dir, ic_file)

        # ROM algo and mode count directory
        if runtype in ["rom", "hyper"]:
            mkdir(rundir)
            rundir = os.path.join(rundir, f"{solve_algo}{nmodes}")
            if runtype == "hyper":
                mkdir(rundir)
                rundir = os.path.join(rundir, sampalgo)
                mkdir(rundir)
                rundir = os.path.join(rundir, gpod_weigher)
                mkdir(rundir)
                finaldir = ""
                if gpod_weigher != "identity":
                    finaldir += f"modes_{nmodes_gpod}_"
                finaldir += f"samp_{sampperc}"
                if seed_qdeim:
                    finaldir += "_qdeim"
                if physrate > 0:
                    finaldir += f"_phys{physrate}"
                if domrate > 0:
                    finaldir += f"_dom{domrate}"
                rundir = os.path.join(rundir, finaldir)

        elif "decomp" in runtype:
            mkdir(rundir)
            # additive/multiplicative directory
            if isadditive:
                rundir = os.path.join(rundir, "additive")
            else:
                rundir = os.path.join(rundir, "multiplicative")
            mkdir(rundir)

            dirname = ""
            for dom_idx in range(ndomains):
                if solve_algo[dom_idx] == "FOM":
                    dirname += "FOM_"
                else:
                    dirname += f"{solve_algo[dom_idx]}{nmodes[dom_idx]}_"
            dirname = dirname[:-1]
            rundir = os.path.join(rundir, dirname)
            if runtype == "hyper_decomp":
                mkdir(rundir)
                rundir = os.path.join(rundir, sampalgo)
                mkdir(rundir)
                rundir = os.path.join(rundir, gpod_weigher)
                mkdir(rundir)
                finaldir = ""
                if gpod_weigher != "identity":
                    assert all([modeval == nmodes_gpod[0] for modeval in nmodes_gpod]), "Can't handle different modes in each hyp-red domain yet"
                    finaldir += f"modes_{nmodes_gpod[0]}_"
                finaldir += f"samp_{sampperc}"
                if seed_qdeim:
                    finaldir += "_qdeim"
                rundir = os.path.join(rundir, finaldir)
                if physrate > 0:
                    rundir += f"_phys{physrate}"
                if domrate > 0:
                    rundir += f"_dom{domrate}"

        if os.path.isdir(rundir) and (not recalc):
            print(f"Not recalculating run at {rundir}")
            continue
        else:
            mkdir(rundir)

        runfile = os.path.join(rundir, "input.yaml")
        with open(runfile, "w") as f:

            f.write(f"loglevel: \"{loglevel}\"\n")
            f.write(f"logtarget: \"{logtarget}\"\n")
            f.write(f"equations: \"{equations}\"\n")
            f.write(f"fluxOrder: {order}\n")
            f.write(f"problemName: \"{problem}\"\n")
            f.write(f"icFlag: {icFlag}\n")
            for param_idx, param in enumerate(params_names_list):
                f.write(f"{param}: {run_list[param_idx]}\n")

            # these may need to be modified?
            f.write(f"meshDirFull: \"{meshdir_full}\"\n")
            f.write(f"odeScheme: \"{scheme}\"\n")
            f.write(f"finalTime: {tf}\n")
            f.write(f"stateSamplingFreq: {sampfreq}\n")
            if "decomp" not in runtype:
                f.write(f"timeStepSize: {dt}\n")
                if ic_index is not None:
                    f.write(f"icFile: \"{ic_file}\"\n")

            if runtype in ["rom", "hyper"]:
                f.write("rom:\n")
                f.write(f"  algorithm: \"{solve_algo}\"\n")
                f.write(f"  numModes: {nmodes}\n")
                f.write(f"  basisFile: \"{basis_root}.bin\"\n")
                f.write(f"  transFile: \"{shift_root}.bin\"\n")

                if runtype == "hyper":
                    f.write("hyper:\n")
                    f.write(f"  meshDirHyper: \"{meshdir_hyper}\"\n")
                    f.write(f"  sampleFile: \"{meshdir_hyper}/sample_mesh_gids.dat\"\n")
                    f.write(f"  stencilFile: \"{meshdir_hyper}/stencil_mesh_gids.dat\"\n")
                    f.write(f"  gpodWeigherType: \"{gpod_weigher}\"\n")
                    if gpod_weigher != "identity":
                        f.write(f"  numModesGpod: {nmodes_gpod}\n")
                        f.write(f"  basisFileGpod: \"{basis_root_gpod}.bin\"\n")

            if "decomp" in runtype:
                f.write("decomp:\n")
                f.write(f"  domainTypes: {solve_algo}\n")
                f.write(f"  timeStepSize: {dt}\n")
                f.write(f"  additive: {isadditive}\n")
                f.write(f"  convStepMax: {conv_step_max}\n")
                if ic_index is not None:
                    f.write(f"  icFileRoot: \"{ic_file}\"\n")

                if runtype in ["rom_decomp", "hyper_decomp"]:
                    f.write(f"  numModes: {nmodes}\n")
                    f.write(f"  basisFileRoot: \"{basis_root}\"\n")
                    f.write(f"  transFileRoot: \"{shift_root}\"\n")
                    if runtype == "hyper_decomp":
                        f.write(f"  sampleFiles: {sampfiles}\n")
                        f.write(f"  gpodWeigherType: \"{gpod_weigher}\"\n")
                        if gpod_weigher != "identity":
                            f.write(f"  gpodBasisRoot: \"{basis_root_gpod}\"\n")
                            f.write(f"  gpodSizeVec: {nmodes_gpod}\n")

        if run:
            print(f"Executing run at {rundir}", flush=True)
            os.chdir(rundir)
            subprocess.call([runner, runfile])
        else:
            print(f"Input file written to {runfile}")

    if not run:
        print("Pass run=True to execute next time")

    # ----- END RUN GENERATION -----

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

    if inputs["runtype"] not in ["fom", "fom_decomp"]:
        solve_algo = inputs["solve_algo"]
    else:
        solve_algo = "FOM"

    # initial conditions
    try:
        ic_index = inputs["ic_index"]
    except KeyError:
        ic_index = None

    # handle ROM inputs
    if "fom" not in inputs["runtype"]:
        nmodes = inputs["nmodes"]
        basis_dir = inputs["basis_dir"]
        basis_file = inputs["basis_file"]
        shift_file = inputs["shift_file"]
    else:
        nmodes = None
        basis_dir = None
        basis_file = None
        shift_file = None

    if "hyper" in inputs["runtype"]:
        sampalgo        = inputs["sampalgo"]
        seed_qdeim      = inputs["seed_qdeim"]
        if (sampalgo != "random") or seed_qdeim:
            nmodes_greedy = inputs["nmodes_greedy"]
        else:
            nmodes_greedy = None
        sampperc        = inputs["sampperc"]

        gpod_weigher    = inputs["gpod_weigher"]
        if gpod_weigher != "identity":
            nmodes_gpod    = inputs["nmodes_gpod"]
            basis_dir_gpod = inputs["basis_dir_gpod"]
        else:
            nmodes_gpod = None
            basis_dir_gpod = None

        if "physrate" in inputs:
            physrate = inputs["physrate"]
        else:
            physrate = 0
        if "domrate" in inputs:
            domrate = inputs["domrate"]
        else:
            domrate = 0
    else:
        sampalgo = None
        nmodes_greedy = None
        sampperc = None
        gpod_weigher = None
        nmodes_gpod = None
        basis_dir_gpod = None
        seed_qdeim = False
        physrate = 0
        domrate = 0

    # handle decomp inputs
    if "decomp" in inputs["runtype"]:
        ndomX = inputs["ndomX"]
        ndomY = inputs["ndomY"]
        overlap = inputs["overlap"]
        isadditive = inputs["isadditive"]
        numprocs = inputs["numprocs"]
        try:
            conv_step_max = inputs["conv_step_max"]
        except KeyError:
            conv_step_max = 10
    else:
        ndomX = None
        ndomY = None
        overlap = None
        isadditive = None
        numprocs = 1
        conv_step_max = None

    gen_runs(
        inputs["equations"],
        inputs["order"],
        inputs["problem"],
        inputs["scheme"],
        inputs["tf"],
        inputs["dt"],
        inputs["sampfreq"],
        inputs["meshroot"],
        inputs["nx"],
        inputs["ny"],
        inputs["runtype"],
        inputs["icFlag"],
        inputs["runner"],
        inputs["loglevel"],
        inputs["logtarget"],
        inputs["outdir_base"],
        phys_params_user=phys_params_user,
        ic_params_user=ic_params_user,
        solve_algo=solve_algo,
        nmodes=nmodes,
        basis_dir=basis_dir,
        basis_file=basis_file,
        shift_file=shift_file,
        ndomX=ndomX,
        ndomY=ndomY,
        overlap=overlap,
        isadditive=isadditive,
        conv_step_max=conv_step_max,
        numprocs=numprocs,
        ic_index=ic_index,
        sampalgo=sampalgo,
        nmodes_greedy=nmodes_greedy,
        sampperc=sampperc,
        gpod_weigher=gpod_weigher,
        nmodes_gpod=nmodes_gpod,
        basis_dir_gpod=basis_dir_gpod,
        seed_qdeim=seed_qdeim,
        physrate=physrate,
        domrate=domrate,
        run=inputs["run"],
    )

