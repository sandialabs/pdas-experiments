# Overview

This repository supplies experiment automation functionality for monolithic and domain-decomposed fluid flow ODEs as constructed by the [pressio-demoapps](https://github.com/Pressio/pressio-demoapps) solver, via the Schwarz decomposition framework [pressio-schwarz](https://github.com/Pressio/pressio-schwarz) and the [pressio](https://github.com/Pressio/pressio) projection-based ROM utilities. This is largely a modification of the [pressio-tutorials](https://github.com/Pressio/pressio-tutorials) methodology, with extensions for the domain decomposition utilities and simplifications where possible.

# Building

Building the runner executable requires (at minimum) a C++17-compatible compiler, and copies of the **pressio**, **pressio-demoapps** (which also includes the **Eigen** linear algebra library), and **pressio-schwarz** source. This repository includes the C++ YAML library, which is built along with the runner executable.

```
export CXX=<path-to-your-CXX-compiler>
export PRESSIO=<path-to-pressio-root>
export PDA=<path-to-pressio-demoapps-root>
export PSCHWARZ=<path-to-pressio-schwarz-root>

git clone git@github.com:sandialabs/pdas-experiments.git
cd pdas-experiments && mkdir build && cd build
cmake -DPRESSIO_SOURCE=${PRESSIO} -DPDA_SOURCE=${PDA} -DPSCHWARZ_SOURCE=${PSCHWARZ} ..
make
```

This will generate two executables: a serial-only version `runner_serial`, and a parallel version (using OpenMP) `runner_omp`.

# Running experiment campaigns

Running experiments is managed through a YAML input file format and the associated experiment execution script found at `python/campaign_runner.py`. To use this utility, you must install the `pschwarz` Python package supplied by **pressio-schwarz** repository. Installing the `pschwarz` package will install most necessary packages, except for the `pyyaml` package which must be installed manually.

To construct a simulation campaign input file, follow the format demonstrated by `cases/jcam-25/2d_euler/case_def.yaml`, including specifying the absolute paths to `runner_serial` and `runner_omp`, the absolute paths to the **pressio-demoapps** (`pdaroot`) and **pressio-schwarz** (`pschwarzroot`) source directories, and the absolute path to the campaign output root directory `caseroot`. Then, given a complete input file `case_def.yaml`, simply execute

```python
python3 python/campaign_runner.py case_def.yaml
```

The campaign outputs will then begin to populate `caseroot`.

# Reproducibility

Input files for experimental campaigns used in publications are stored in `cases/`. Citation information for the relevant publication is noted in the README of the subdirectories. Results may slightly change due to modifications of **pressio** and **pressio-demoapps**. Additionally, Python scripts for generating the images presented in the publications can be found in the `figures/` directory of each case subdirectory.