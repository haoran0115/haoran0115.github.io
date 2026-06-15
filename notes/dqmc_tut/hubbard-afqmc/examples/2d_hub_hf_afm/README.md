# 2D half-filled Hubbard AFM example

This example runs the minimal projector DQMC code for a $6 \times 6$ half-filled Hubbard model with

- $U = 4.0$
- `ham_mu = 0.0`
- `hs_channel = :SU2`
- `Nbin = 20`
- `Nsweep = 100`

The Julia source files are symlinked to the main codebase. Only `parameters.jl` and `plot.py` are local to this example.

The current projector workflow is canonical: `Npart = Ns / 2` is fixed internally, so `ham_mu` enters the one-body Hamiltonian but does not tune the filling in this example.

## Run

From this directory:

```bash
julia main.jl
python plot.py
```

Optional MPI run:

```bash
mpiexecjl -np 4 julia main.jl
python plot.py
```

## Output

`julia main.jl` writes raw bin data under `data/`, runs the built-in analysis path, and writes analyzed observables under `analysis/`.

`python plot.py` reads `analysis/data.h5`, prints a short summary of the AFM structure factor, and saves

- `afm_structure_factor_heatmap.png`
