# Projector DQMC code for the Hubbard model

This directory contains a small Julia implementation of projector DQMC for the Hubbard model.

The current workflow is:

- edit [parameters.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/parameters.jl)
- run [main.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/main.jl)
- inspect raw data under `data/` and analyzed data under `analysis/`

## Main files

- [main.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/main.jl): runs the Monte Carlo simulation
- [hubbard-afqmc.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/hubbard-afqmc.jl): core data structures, updates, measurements, and analysis helpers
- [parameters.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/parameters.jl): model and QMC parameters
- [analysis.jl](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/analysis.jl): reruns the built-in postprocessing on saved data

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

- [examples/2d_hub_hf_afm](/Users/shiroha/Documents/websites/haoran0115.github.io/notes/dqmc_tut/hubbard-afqmc/examples/2d_hub_hf_afm): runs a $6 \times 6$ half-filled Hubbard model and plots the AFM structure factor in momentum space.
