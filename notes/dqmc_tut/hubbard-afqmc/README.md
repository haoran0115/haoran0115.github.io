# Projector DQMC code for the Hubbard model

This directory contains a small Julia implementation of projector DQMC for the Hubbard model.

The current workflow is:

- edit [parameters.jl](parameters.jl)
- run [main.jl](main.jl)
- inspect raw data under `data/` and analyzed data under `analysis/`

## Main files

- [main.jl](main.jl): runs the Monte Carlo simulation
- [hubbard-afqmc.jl](hubbard-afqmc.jl): core data structures, updates, measurements, and analysis helpers
- [parameters.jl](parameters.jl): model and QMC parameters
- [analysis.jl](analysis.jl): reruns the built-in postprocessing on saved data

## Run

From this directory:

```bash
julia main.jl
```

Optional MPI run:

```bash
mpiexecjl -np 4 julia main.jl
```

If you already have `data/data.h5` and want to rerun the postprocessing only:

```bash
julia analysis.jl
```

## Examples

For real runnable examples that use this code, see:

- [examples/2d_hub_hf_afm](examples/2d_hub_hf_afm): runs a $6 \times 6$ half-filled Hubbard model and plots the AFM structure factor in momentum space.
