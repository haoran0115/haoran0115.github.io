---
layout: default
title: DQMC Tutorial
---

# An Introduction to Determinant Quantum Monte Carlo for Interacting Fermion Systems

## Scope
Determinant quantum Monte Carlo (DQMC), often called auxiliary-field quantum Monte Carlo in the lattice-fermion context, is a numerically exact method for computing thermal ($T > 0$) or ground state ($T = 0$) properties of interacting fermion systems. The central idea is to avoid explicit diagonlization over the exponentially large many-body Hilbert space, by rewriting the interacting quantum problem as a stochastic summation over auxiliary-field configurations.

Two senarios are discussed here
1. **Finite-temperature DQMC**, where one want to compute finte-temperature properties
   $$
   \braket{\hat O}
   = \frac{\text{Tr}\, \hat O e^{-\beta \hat H}}
   {\text{Tr}\, e^{-\beta \hat H}}
   $$
2. **Projector DQMC**, where one starts from a trial Slater determinant $\ket{\Psi_T}$ with $N_\text{part}$ number of particles and computes
   $$
   \braket{\hat {O}}
   = \lim_{\Theta\to\infty}
   \frac{\braket{\Psi_T|e^{-\Theta \hat H}\hat O e^{-\Theta \hat H}|\Psi_T}}
   {\braket{\Psi_T|e^{-2\Theta \hat H}|\Psi_T}}
   $$


## Many-body numerics
Before introducing determinant quantum Monte Carlo, it is useful to clarify why generic interacting fermion problems are difficult, and to introduce the notation used in these notes.

### Hamiltonian of interacting fermions
A general number-conserving two-body fermion can be written as
$$
    \hat H = \hat H_0 + \hat H_I,\quad
    H_0 = \sum_{ij} T_{ij}c_i^\dagger c_j,\quad
    H_I = \frac{1}{2} \sum_{ijkl} V_{ijkl} c_i^\dagger c_j^\dagger c_l c_k
$$
where $\hat H_0$ is the non-interacting part and $\hat H_I$ is the interacting part. $i,j,k,l = 1,2,\dots,N$ labels different degrees of freedom in the system, e.g. for Hubbard model on lattices, $(i,\sigma)$ labels spatial and spin degrees of freedom.
$$
    \hat H = -t \sum_{\braket{ij},\sigma} (c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.}) + U\sum_i n_{i\uparrow} n_{i\downarrow}
$$

If only $H_0$ is non-zero, one can directly diagonalize the one-body matrix $T = U^\dagger \lambda U$ to get well-defined energy levels. The systems are just well-defined single particles, and the computaitonal complexity to diagonalize this $N\times N$ matrix is $O(N^3)$.
$$
    H_0 = \sum_{ij} T_{ij}c_i^\dagger c_j = \sum_{ijk} (c_i U_{ik}^\dagger) \lambda_i  (U_{kj} c_j^\dagger) = \sum_k \lambda_k d_k^\dagger d_k
$$

Howerver if the interaction $H_I$ is non-zero, there are no well-defined single-particle energy levels and one must explicitly write down all basis to explicitly diagonlize the Hamiltonian.
The number of basis (i.e. the size of the Hilbert space) scales exponentially and it makes impossible to solve large system using exact diagonalization.

For example, a 4-site spinless fermion system have $2^4 = 16$ basis vectors
$$
    \mathcal{B} = \{
        \ket{0000},\ket{0001}, \ket{0010}, ..., \ket{1110}, \ket{1111}
    \}
$$
where $0/1$ means zero/one fermion occupation.

This exponential growth of the many-body Hilbert space is the basic bottleneck. DQMC avoids explicitly enumerating and diagonalizing the many-body Hamiltonian in this basis. Instead, it rewrites the interacting fermion problem as a stochastic sum over auxiliary-field configurations, where each configuration corresponds to a non-interacting fermion problem.

## From MCMC to DQMC
This section begins with a brief review of Markov-chain Monte Carlo (MCMC) sampling in the classical Ising model, and then uses this example to motivate the stochastic formulation underlying determinant quantum Monte Carlo.

### MCMC for Ising model
For the classical Ising model on a $N_s$-sites lattice
$$
H_\text{Ising} = -J\sum_{\braket{ij}} s_i s_j,
\quad s_i=\pm 1,
$$
the partition function and thermal expectation value of an observable $O[s]$ can be expressed by a summation over all possible spin configurations $\mathcal{C} = \{s_i\}$
$$
Z=\sum_{\{s_i\}} e^{-\beta H_{\mathrm{Ising}}[\{s_i\}]} \equiv \sum_\mathcal{C} W[\mathcal{C}],\quad
\langle O\rangle
= \frac{1}{Z}\sum_\mathcal{C} O[\mathcal{C}] W[\mathcal{C}]
$$
Instead of summing over all $2^{N_s}$ configurations, we sample configurations with probability (we can do this since $W[\mathcal{C}]$ is always nonnegative!)
$$
P[\mathcal{C}] =  \frac{W[\mathcal{C}]}{Z}
$$

A simple Metropolis update proposes a local spin flip $s_i\to s_i' = -s_i$, computes the weight ratio and acceptance probability
$$
R = \frac{W[\mathcal{C}']}{W[\mathcal{C}]} = e^{-\beta(H[\{s'_i\}]-H[\{s_i\}])},\quad A(s_i\to s_i')=\min(1,R)
$$
The result is a Markov chain whose long-time distribution is the desired Boltzmann distribution, provided the update is ergodic and satisfies detailed balance.

Here we reformulate the MCMC algorithm above using a more general language
1. From current existing configuration $\mathcal{C}$, propose new configuration $\mathcal{C}'$ (e.g. through local spin flip)
1. Calculate the acceptance probablity $R$ and accept or reject 
1. This process effectively samples $\mathcal{C}$ from a probability distribution $P[\mathcal{C}] = \frac{W[\mathcal{C}]}{Z}$

The Ising model only contains nearest neighbor couplings, so the cost of evaluating a single Monte Carlo acceptance probability is $O(1)$. A single sweep (which enumerates all $N$ spin flips) will give $O(N)$ complexity.

The DQMC problem will have the same mathematical sturcture
$$
\langle \hat O\rangle
=\frac{\sum_s W[\mathcal{C}] \braket{\hat O}_\mathcal{C}}{\sum_s W[\mathcal{C}]}
$$
but the configuration $s$ will be a space-time auxiliary field, and the weight $W[\mathcal{C}]$ will be a fermion determinant.

### DQMC targets
Here we reiterate what we want to compute.
For finite temperature
$$
Z_\beta = \text{Tr}\,\,e^{-\beta \hat H},\quad
\langle \hat O\rangle
=\frac{\text{Tr}\,\,\hat O e^{-\beta\hat H}}
{Z_\beta}
$$
For ground-state projection
$$
Z_\Theta=\braket{\Psi_T|e^{-2\Theta\hat H}|\Psi_T},\quad
\langle \hat O\rangle_0
=\lim_{\Theta\to\infty}
\frac{\braket{\Psi_T|e^{-\Theta\hat H}\hat O e^{-\Theta\hat H}|\Psi_T}}
{Z_\Theta}
$$
Here $\ket{\Psi_T}$ is a trial state with nonzero ground-state overlap,
$$
\braket{\Psi_T|\Psi_0}\ne0
$$

Recall
* if $\hat H$ is non-interacting, evaluate these quantites will only takes up to $O(N^3)$ scaling
* if $\hat H$ is interacting, since the Hilbert space exponentially scales with system size, it would be exponentially hard to compute $ e^{-\beta\hat H}$ and $e^{-\Theta\hat H}$

Both finite-temperature and projector DQMC prevent this use the same next step, which is the central ider of QMC
> Discretize the interacting propagator $e^{-\beta\hat H}$ and $e^{-2\Theta\hat H}$ into stochastic sums of non-interacting propagators

### Stochastic representation of propagator
All stories of DQMC starts from the following identity
$$
    e^{A^2 / 2} = \frac{1}{\sqrt{2\pi}}\int dq~ e^{-(q-A)^2 / 2 + A^2 / 2} 
    = \frac{1}{\sqrt{2\pi}} \int dq~ e^{-q^2 / 2}  e^{qA}
    = \mathbb{E}[e^{qA}]
$$
$q$ is the so-called auxiliary-field: it decouples $e^{A^2}$ into a expecation under stochastic variable $q$ which follows a $\mathcal{N}(0,1)$ Gaussian distribution.
This equation also holds when $A$ is an operator.
While in practice, one usually approximate this integral use a four-components field $l = \pm 1, \pm 2$ 
$$
    e^{A^2 / 2} = \frac{1}{\sqrt{2\pi}} \int dq~ e^{-q^2 / 2}  e^{qA}
    \approx \sum_{l=\pm 1, \pm 2} \gamma_l e^{\eta_l A}
$$
(might give a picture here)
Values of $\gamma_l,\eta_l$ can be found from Gauss–Hermite quadrature
$$
    \gamma_{\pm 1} = \frac{1}{4}\left( 1 - \frac{\sqrt{6}}{3} \right),\quad
    \gamma_{\pm 2} = \frac{1}{4}\left( 1 + \frac{\sqrt{6}}{3} \right)\\
    \eta_{\pm 1} = \pm \sqrt{2(3 + \sqrt{6})},\quad
    \eta_{\pm 2} = \pm \sqrt{2(3 - \sqrt{6})}
$$
$q/l$ is called continous/discrete auxiliary field.

Suppose, the interaction term $\hat H_I = \sum_i \lambda_i \hat A_i^2$ can be written in terms of a summation of squares of non-interacting operator $\hat A_i$, its propagator with propagation length $\Delta\tau$ can also be represented using the stochastic summation (explain what is a propagator)
$$
    e^{-\Delta\tau \hat H_I} = \prod_i e^{-\lambda\Delta\tau \hat A_i^2}
    = \prod_i \left[
        \sum_{l_i = \pm1, \pm2}
        \gamma_{l_i} e^{\eta_{l_i} \sqrt{-\lambda_i\Delta\tau} \hat A_i}
    \right] + O(\Delta\tau^4)
$$
where each term in the summation is an non-interacting propagator!
For example, the Hubbard interaction can be re-write in the following two ways (up to a chemical potential shift)
$$
    \hat H_I = \pm \frac{U}{2}\sum_i 
    (c_{i\uparrow}^\dagger c_{i\uparrow} \pm 
    c_{i\downarrow}^\dagger c_{i\downarrow})^2,\quad
    \hat A_i = c_{i\uparrow}^\dagger c_{i\uparrow} \pm 
    c_{i\downarrow}^\dagger c_{i\downarrow}
$$

Suzuki-Trotter decomposition enables us to represent a short and interacting propagator $e^{-\Delta\tau\hat H}$ 
$$
    e^{-\Delta\tau\hat H} = e^{-\Delta\tau\hat H_0 - \Delta\tau\hat H_I}
    = e^{-\Delta\tau\hat H_0} e^{-\Delta\tau\hat H_I} + O(\Delta\tau)
$$

Then, we can assemble a long propagator with total propagation length $N_\tau \Delta\tau = \beta$ or $2\Theta$
$$
\begin{aligned}
    e^{-\beta\hat H} \text{ or } e^{-2\Theta\hat H} 
    &= 
    \underbrace{
        e^{-\Delta\tau \hat H}e^{-\Delta\tau \hat H}\cdots
        e^{-\Delta\tau \hat H}
    }_{N_\tau\text{ copies in total}}\\
    &= \prod_{t = 1}^{N_\tau} \left\{
        e^{-\Delta\tau\hat H_0}
        \left[
            \sum_{l_{t,i} = \pm1, \pm2}
            \gamma_{l_i} e^{\eta_{l_i} \sqrt{-\lambda_i\Delta\tau} \hat A_i}
        \right]
    \right\} + O(\Delta\tau)\\
    &= \sum_{\{l_{t,i}\}} \prod_{t = 1}^{N_\tau} \left\{
        e^{-\Delta\tau\hat H_0}
        \left[
            \gamma_{l_i} e^{\eta_{l_i} \sqrt{-\lambda_i\Delta\tau} \hat A_i}
        \right]
    \right\}  + O(\Delta\tau)\\
    &= \sum_{\mathcal{C}} U_\mathcal{C} + O(\Delta\tau)
\end{aligned}
$$
where in the expression above
1. $t = 1,2,\dots,N_\tau$ labels imaginary time discretization index
1. $i$ labels interaction vertices $A_i$, e.g. for Hubbard model, each lattice site have one interaction vertex.
1. Auxiliary field $l_{t,i}$ is also labeled by these two indices
1. We call an auxilary field configuration $\mathcal{C} = \{l_{t,i}\}$, and called the corresponding propagator $U_\mathcal{C}$.

*Remark.* It could be more natural working under path integral formalism
$$
    \int\mathcal{D}[\bar c, c]
    ~e^{\int d\tau~ -\bar c (\partial_t - T) c + A^2}
    = C \int\mathcal{D}[\bar c, c] \mathcal{D}[q]
    ~e^{\int d\tau~ -\bar c (\partial_t - T) c + qA - q^2/2}
$$
$q = q_i(\tau)$ is the auxiliary field. $q = \braket{A}$ is the so-called mean-field/saddle-point approximation.

In summary, the stochastic representation of propagator is achieved by two steps
1. Suzuki-Trotter decomposition of a long propagator
2. Discrete Hubbard-Stratonovich transformation which convert interacting propagators into summations of non-interacting propagators

### Fermion weights
The partition function can then be formulated by a summation over auxiliary field configurations
$$
\begin{aligned}
    Z_\beta &= \sum_\mathcal{C} W[\mathcal{C}],\quad W[\mathcal{C}] = \text{Tr} U_\mathcal{C} \\
    Z_\Theta &= \sum_\mathcal{C} W[\mathcal{C}],\quad W[\mathcal{C}] = \braket{\Psi_T| U_\mathcal{C} |\Psi_T}
\end{aligned}
$$
This equation looks almost the same as Ising model, except we don't konw if the weight $W[\mathcal{C}]$ is nonnegative or not. If $W[\mathcal{C}]$ is nonnegative, we can build an exactly same MCMC algorithm comparing to 

Since $U_\mathcal{C}$ represents a squence of $N_\tau$ non-interacting propagators, the complexity for computing
* $W[\mathcal{C}]$ is $O(N_\tau N^3)$

It can also be shown that the comutational complexity for a single flipping sweep can be reduced to the same.

Here, we summarize the difference between the Ising (classical) MC and DQMC
|    | Ising MC | DQMC |
|:-- | :-- | :-- |
| Field config $\mathcal{C}$ | $\{s_i\}$: space | $\{l_{t,i}\}$: time, space |
| Update complexity | $O(N)$ | $O(N_\tau N^3)$ | 
| Weight | $W[\mathcal{C}]\ge 0$ | $W[\mathcal{C}]\ngeq 0$ |

## Fermion sign problem
The determinant weight $W[\mathcal{C}]$ is not guaranteed to be positive. It can be negative or complex. At this stage $P[\mathcal{C}] = W[\mathcal{C}]/Z$ is no longer a well-defined probability. Instead, we do the Monte Carlo for the following well-defined probabilistic distribution
$$
\begin{aligned}
    Z' &= \sum_\mathcal{C} |\mathrm{Re}\,W[\mathcal{C}]|, \\
    P[\mathcal{C}] &= \frac{|\mathrm{Re}\,W[\mathcal{C}]|}{Z'}, \\
    \text{sign}[\mathcal{C}] &= \frac{W[\mathcal{C}]}{|\mathrm{Re}\,W[\mathcal{C}]|}, \\
    \braket{\hat O}
    &= \frac{\sum_\mathcal{C} |\mathrm{Re}\,W[\mathcal{C}]| \, \text{sign}[\mathcal{C}] \, \braket{\hat O}_\mathcal{C}}
            {\sum_\mathcal{C} |\mathrm{Re}\,W[\mathcal{C}]| \, \text{sign}[\mathcal{C}]}
     = \frac{\braket{\hat O \, \text{sign}}_{|\mathrm{Re}\,W|}}
            {\braket{\text{sign}}_{|\mathrm{Re}\,W|}}
\end{aligned}
$$
The MCMC algorithm is then
1. From current existing configuration $\mathcal{C}$, propose new configuration $\mathcal{C}'$ (e.g. through local spin flip)
1. Calculate the acceptance probability $R$ and accept or reject 
    $$
        R = \frac{|\mathrm{Re}\,W[\mathcal{C'}]|}{|\mathrm{Re}\,W[\mathcal{C}]|}
    $$
1. This process effectively samples $\mathcal{C}$ from a probability distribution $P[\mathcal{C}] = \frac{|\mathrm{Re}\,W[\mathcal{C}]|}{Z'}$

There are also other reweighting schemes such as taking $\mathrm{Re} W[\mathcal{C}]$ rather than absolute value.

## Pratical implementation
TBD

## Minimal Implementation
Source code in Julia are under [hubbard-afqmc](https://github.com/haoran0115/haoran0115.github.io/tree/main/hubbard-afqmc) folder
* `hubbard-afqmc.jl`: utility functions
* `main.jl`: main QMC loop
* `parameters.jl`: simulation paramters
* `seeds.jl`: random seed

How to run the code: modify the simulation parameters in `parameters.jl` and run with
```bash
# with single markov chain
julia main.jl

# with N markov chains, e.g. N = 4 for the following command
mpiexecjl -np 4 julia main.jl
```

Simulation results will be written under `data/` and analyzed results will be written under `analysis/`.

## Reference and useful resources

Tutorials
* [TowardQMC](https://www.youtube.com/playlist?list=PLheYERt_Ks3wHU_Fa7MVGH5I31pJf4SLH): video introduction for finite-$T$ formalism
* [DQMC](https://quantummc.xyz/teaching/dqmc/): concise introduction for finite-$T$ DQMC using Hubbard model as an example, written Chinese

Technical walk-through
* [World line and determinantal Quantum Monte
Carlo methods for spins, phonons, and
electrons](https://pawn.physik.uni-wuerzburg.de/~assaad/Reprints/assaad_evertz.pdf): part of this tutorial and most of the minimal implementation is mainly based on this tutorial

Software packages
* [ALF](https://git.physik.uni-wuerzburg.de/ALF/ALF): very comprehensive DQMC/AFQMC package in Fortran, with complete [documentation](documentation), minimal tutorials, and [workshop recordings&tutorials](https://git.physik.uni-wuerzburg.de/ALF/ALF_Tutorial/-/tree/master/Presentations?ref_type=heads)
* [SmoQyDQMC](https://smoqysuite.github.io/SmoQyDQMC.jl/stable/): HMC (an determinant-free QMC) for electron-electron and electron-phonon problems
