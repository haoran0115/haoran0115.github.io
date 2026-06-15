---
layout: default
title: DQMC Tutorial
---

# An Introduction to Determinant Quantum Monte Carlo for Interacting Fermion Systems
{:.no_toc}

* TOC
{:toc}

## Scope
Determinant quantum Monte Carlo (DQMC), often called auxiliary-field quantum Monte Carlo in the lattice-fermion context, is a numerically exact method for computing thermal ($T > 0$) or ground state ($T = 0$) properties of interacting fermion systems. The central idea is to avoid explicit diagonalization over the exponentially large many-body Hilbert space, by rewriting the interacting quantum problem as a stochastic summation over auxiliary-field configurations.

Two scenarios are discussed here:
1. **Finite-temperature DQMC**, where one wants to compute finite-temperature properties

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

A slide version is available [here](/notes/dqmc_tut/slides/main.pdf).


## Many-body numerics
Before introducing determinant quantum Monte Carlo, it is useful to clarify why generic interacting fermion problems are difficult, and to introduce the notation used in these notes.

A general number-conserving two-body fermion Hamiltonian can be written as

$$
    \hat H = \hat H_0 + \hat H_I,\quad
    \hat H_0 = \sum_{ij} T_{ij}c_i^\dagger c_j,\quad
    \hat H_I = \frac{1}{2} \sum_{ijkl} V_{ijkl} c_i^\dagger c_j^\dagger c_l c_k
$$

where $\hat H_0$ is the non-interacting part and $\hat H_I$ is the interacting part. $i,j,k,l = 1,2,\dots,N$ labels different degrees of freedom in the system, e.g. for Hubbard model on lattices, $(i,\sigma)$ labels spatial and spin degrees of freedom.

$$
    \hat H = -t \sum_{\braket{ij},\sigma} (c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.}) + U\sum_i n_{i\uparrow} n_{i\downarrow}
$$

If only $\hat H_0$ is non-zero, one can directly diagonalize the one-body matrix $T = U \lambda U^\dagger$ to obtain well-defined single-particle energy levels. The computational cost of diagonalizing this $N\times N$ matrix is $O(N^3)$.

$$
    \hat H_0 = \sum_{ij} T_{ij}c_i^\dagger c_j = \sum_{ijk} (c_i^\dagger U_{ik}) \lambda_k  (U_{kj}^\dagger c_j) = \sum_k \lambda_k d_k^\dagger d_k
$$

However, if the interaction $\hat H_I$ is non-zero, there are no well-defined single-particle energy levels and one must work in the full many-body Hilbert space. Explicitly diagonalizing the Hamiltonian requires writing down all basis states.
The number of basis states (i.e. the size of the Hilbert space) scales exponentially, making it impossible to solve large systems using exact diagonalization.

For example, a 4-site spinless fermion system has $2^4 = 16$ basis vectors

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

the partition function and thermal expectation value of an observable $O[s]$ can be expressed by a summation over all possible spin configurations $\mathcal{C} = \lbrace s_i \rbrace$

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
1. From the current configuration $\mathcal{C}$, propose a new configuration $\mathcal{C}'$ (e.g. through a local spin flip).
1. Calculate the weight ratio $R$ and accept or reject with probability $\min(1, \lvert \mathrm{Re}\,R \rvert)$.
1. This process effectively samples $\mathcal{C}$ from a probability distribution $P[\mathcal{C}] = \frac{W[\mathcal{C}]}{Z}$

The Ising model only contains nearest neighbor couplings, so the cost of evaluating a single Monte Carlo acceptance probability is $O(1)$. A single sweep (which enumerates all $N$ spin flips) will give $O(N)$ complexity.

The DQMC problem will have the same mathematical structure

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

Recall:
* If $\hat H$ is non-interacting, evaluating these quantities takes only $O(N^3)$ effort.
* If $\hat H$ is interacting, the exponential growth of the Hilbert space makes computing $e^{-\beta\hat H}$ or $e^{-\Theta\hat H}$ exponentially hard.

Both finite-temperature and projector DQMC circumvent this bottleneck through the same central idea:
> Discretize the interacting propagator $e^{-\beta\hat H}$ and $e^{-2\Theta\hat H}$ into stochastic sums of non-interacting propagators

### Stochastic representation of propagator
DQMC begins with the following identity

$$
    e^{A^2 / 2} = \frac{1}{\sqrt{2\pi}}\int dq~ e^{-(q-A)^2 / 2 + A^2 / 2} 
    = \frac{1}{\sqrt{2\pi}} \int dq~ e^{-q^2 / 2}  e^{qA}
    = \mathbb{E}[e^{qA}]
$$

$q$ is the so-called auxiliary-field: it decouples $e^{A^2}$ into an expectation under a stochastic variable $q$ which follows a $\mathcal{N}(0,1)$ Gaussian distribution.
This equation also holds when $A$ is an operator.
In practice one usually approximates this integral with a discrete four-component field $l = \pm 1, \pm 2$:

$$
    e^{A^2 / 2} = \frac{1}{\sqrt{2\pi}} \int dq~ e^{-q^2 / 2}  e^{qA}
    \approx \sum_{l=\pm 1, \pm 2} \gamma_l e^{\eta_l A}
$$

Values of $\gamma_l,\eta_l$ can be found from Gauss–Hermite quadrature

$$
\begin{aligned}
    \gamma_{\pm 1} = \frac{1}{4}\left( 1 - \frac{\sqrt{6}}{3} \right),\quad
    \gamma_{\pm 2} = \frac{1}{4}\left( 1 + \frac{\sqrt{6}}{3} \right)\\
    \eta_{\pm 1} = \pm \sqrt{2(3 + \sqrt{6})},\quad
    \eta_{\pm 2} = \pm \sqrt{2(3 - \sqrt{6})}
\end{aligned}
$$

The variable $q$ (continuous) or $l$ (discrete) is called the auxiliary field.

Now suppose the interaction can be written as a sum of squares, $\hat H_I = \sum_i \lambda_i \hat A_i^2$, where each $\hat A_i = \sum_{ab} [A_i]_{ab}\, c_a^\dagger c_b + a_i$ is a Hermitian fermion bilinear plus a c-number constant $a_i$, and $\lambda_i \in \mathbb{R}$. Its short-time propagator can then be represented using the stochastic summation:
$$
    e^{-\Delta\tau \hat H_I} = \prod_i e^{-\lambda_i\Delta\tau \hat A_i^2}
    = \prod_i \left[
        \sum_{l_i = \pm1, \pm2}
        \gamma_{l_i} e^{\eta_{l_i} \sqrt{-\lambda_i\Delta\tau}\, \hat A_i}
    \right] + O(\Delta\tau^4)
$$

where each term in the summation is a non-interacting propagator. Note that if $\hat A_i$ contains a constant $a_i$, it contributes a scalar factor $e^{\eta_l \sqrt{-\lambda_i\Delta\tau}\,a_i}$ to the weight.
For example, the Hubbard interaction admits two standard decouplings (up to a chemical potential shift):

- **Spin channel:** $U n_{i\uparrow} n_{i\downarrow} = -\frac{U}{2}(n_{i\uparrow} - n_{i\downarrow})^2 + \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow})$, so $\lambda_i = -U/2$ and $\hat A_i = n_{i\uparrow} - n_{i\downarrow}$ with $a_i = 0$.
- **Charge ($SU(2)$) channel:** $U n_{i\uparrow} n_{i\downarrow} = \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow} - 1)^2 - \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow}) + \frac{U}{2}$, so $\lambda_i = U/2$ and $\hat A_i = n_{i\uparrow} + n_{i\downarrow} - 1$ with $a_i = -1$. The last two terms are a chemical potential shift and an irrelevant constant.

The Suzuki–Trotter decomposition splits the short-time propagator as

$$
    e^{-\Delta\tau\hat H} = e^{-\Delta\tau\hat H_0 - \Delta\tau\hat H_I}
    = e^{-\Delta\tau\hat H_0} e^{-\Delta\tau\hat H_I} + O(\Delta\tau^2).
$$

Assembling $N_\tau$ such factors, the total Trotter error accumulates to $O(N_\tau \Delta\tau^2) = O(\Delta\tau)$.

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
            \gamma_{l_{t,i}} e^{\eta_{l_{t,i}} \sqrt{-\lambda_i\Delta\tau}\, \hat A_i}
        \right]
    \right\} + O(\Delta\tau)\\
    &= \sum_{\{l_{t,i}\}} \prod_{t = 1}^{N_\tau} \left\{
        e^{-\Delta\tau\hat H_0}
        \left[
            \gamma_{l_{t,i}} e^{\eta_{l_{t,i}} \sqrt{-\lambda_i\Delta\tau}\, \hat A_i}
        \right]
    \right\}  + O(\Delta\tau)\\
    &= \sum_{\mathcal{C}} U_\mathcal{C} + O(\Delta\tau)
\end{aligned}
$$

where in the expression above
1. $t = 1,2,\dots,N_\tau$ labels imaginary time discretization index
1. $i$ labels interaction vertices, e.g. for Hubbard model, each lattice site has one interaction vertex.
1. Auxiliary field $l_{t,i}$ is also labeled by these two indices
1. We call $\mathcal{C} = \lbrace l_{t,i} \rbrace$ an auxiliary-field configuration, and denote the corresponding propagator by $U_\mathcal{C}$.

*Remark.* It is more natural to work within the path integral formalism

$$
    \int\mathcal{D}[\bar c, c]
    ~e^{\int d\tau~ -\bar c (\partial_t - T) c + A^2}
    = C \int\mathcal{D}[\bar c, c] \mathcal{D}[q]
    ~e^{\int d\tau~ -\bar c (\partial_t - T) c + qA - q^2/2}
$$

Here $q_i(\tau)$ is the auxiliary field; setting $q = \braket{A}$ recovers the mean-field (saddle-point) approximation.

In summary, the stochastic representation of the propagator is achieved by two steps:
1. Suzuki–Trotter decomposition of the long propagator into $N_\tau$ short-time factors.
2. Discrete Hubbard–Stratonovich transformation, which converts each interacting short-time factor into a sum over non-interacting propagators.

### Fermion weights
The partition function can then be formulated by a summation over auxiliary field configurations

$$
\begin{aligned}
    Z_\beta &= \sum_\mathcal{C} W[\mathcal{C}],\quad W[\mathcal{C}] = \text{Tr}\, U_\mathcal{C} \\
    Z_\Theta &= \sum_\mathcal{C} W[\mathcal{C}],\quad W[\mathcal{C}] = \braket{\Psi_T| U_\mathcal{C} |\Psi_T}
\end{aligned}
$$

At the single-particle level, we define the post-HS vertex operator $\hat V_i = \sqrt{-\lambda_i\Delta\tau}\, \hat A_i$ and its matrix $V_i = \sqrt{-\lambda_i\Delta\tau}\, A_i$. The constant in $\hat V_i$ is $v_i = \sqrt{-\lambda_i\Delta\tau}\, a_i$, giving a scalar factor $e^{\eta_{l_{t,i}} v_i}$ in each short-time factor $U_t$. For one imaginary-time slice $t$, the single-particle matrix is

$$
    B_t
    = e^{-\Delta\tau T}
    \prod_i e^{\eta_{l_{t,i}} V_i},
$$

and the scalar prefactor accumulated from all slices is $\prod_{t,i} \gamma_{l_{t,i}} e^{\eta_{l_{t,i}} v_i}$. For the full auxiliary-field history $\mathcal{C} = \lbrace l_{t,i} \rbrace$ define

$$
    B_\mathcal{C}
    = B_{N_\tau} B_{N_\tau-1} \cdots B_1.
$$

In particular, the finite-temperature weight can be written as

$$
    W[\mathcal{C}]
    = \text{Tr}\, U_\mathcal{C}
    = \left(\prod_{t,i}\gamma_{l_{t,i}} e^{\eta_{l_{t,i}} v_i}\right)
      \det[I + B_\mathcal{C}].
$$

If the trial state $\ket{\Psi_T} = \prod_j \left(\sum_i c_i^\dagger P_{ij} \right) \ket{0}$ is a Slater determinant represented by an $N\times N_\text{part}$ orbital matrix $P$, then in projector DQMC the corresponding weight is

$$
    W[\mathcal{C}]
    = \braket{\Psi_T| U_\mathcal{C} |\Psi_T}
    = \left(\prod_{t,i}\gamma_{l_{t,i}} e^{\eta_{l_{t,i}} v_i}\right)
      \det[P^\dagger B_\mathcal{C} P].
$$

This equation looks almost the same as the Ising model, except we do not yet know whether the weight $W[\mathcal{C}]$ is nonnegative. If $W[\mathcal{C}]$ is nonnegative, we can apply the same MCMC algorithm as in the Ising model.

Since $U_\mathcal{C}$ represents a sequence of $N_\tau$ non-interacting propagators, the complexity for computing $W[\mathcal{C}]$ is $O(N_\tau N^3)$.

It can also be shown that the computational complexity for a single flipping sweep can be reduced to the same order.

Here we summarize the differences between classical Ising MC and DQMC:

|    | Ising MC | DQMC |
|:-- | :-- | :-- |
| Field config $\mathcal{C}$ | $\lbrace s_i \rbrace$: space | $\lbrace l_{t,i} \rbrace$: time, space |
| Update complexity | $O(N)$ | $O(N_\tau N^3)$ | 
| Weight | $W[\mathcal{C}]\ge 0$ | $W[\mathcal{C}]\ngeq 0$ |

## Fermion sign problem
The determinant weight $W[\mathcal{C}]$ is not guaranteed to be positive. It can be negative or complex. At this stage $P[\mathcal{C}] = W[\mathcal{C}]/Z$ is no longer a well-defined probability. Instead, we do the Monte Carlo for the following well-defined probabilistic distribution

$$
\begin{aligned}
    Z' &= \sum_\mathcal{C} \lvert \mathrm{Re}\,W[\mathcal{C}] \rvert, \\
    P[\mathcal{C}] &= \frac{\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert}{Z'}, \\
    \text{sign}[\mathcal{C}] &= \frac{W[\mathcal{C}]}{\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert}, \\
    \braket{\hat O}
    &= \frac{\sum_\mathcal{C} \lvert \mathrm{Re}\,W[\mathcal{C}] \rvert \, \text{sign}[\mathcal{C}] \, \braket{\hat O}_\mathcal{C}}
            {\sum_\mathcal{C} \lvert \mathrm{Re}\,W[\mathcal{C}] \rvert \, \text{sign}[\mathcal{C}]}
     = \frac{\braket{\hat O \, \text{sign}}_{\lvert \mathrm{Re}\,W \rvert}}
            {\braket{\text{sign}}_{\lvert \mathrm{Re}\,W \rvert}}
\end{aligned}
$$

The MCMC algorithm is then
1. From the current configuration $\mathcal{C}$, propose a new configuration $\mathcal{C}'$ (e.g. through a local spin flip).
1. Calculate the weight ratio $R$ and accept or reject with probability $\min(1, \lvert \mathrm{Re}\,R \rvert)$.

    $$
        R = \frac{\lvert \mathrm{Re}\,W[\mathcal{C'}] \rvert}{\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert}
    $$

1. This process effectively samples $\mathcal{C}$ from a probability distribution
$P[\mathcal{C}] = \frac{\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert}{Z'}$

## Practical implementation of projector DQMC

This section walks through the practical structure of a projector DQMC simulation for the Hubbard model. The presentation follows the Julia implementation in [`hubbard-afqmc`](https://github.com/haoran0115/haoran0115.github.io/tree/main/notes/dqmc_tut/hubbard-afqmc).

### Hamiltonian

We work with the Hubbard Hamiltonian

$$
\hat H = \hat H_0 + \hat H_I, \quad
\hat H_0 = -t\sum_{\langle ij\rangle,\sigma} c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.},
\quad
\hat H_I = U\sum_i n_{i\uparrow} n_{i\downarrow}
$$

$T_{ij}$ is built from $\hat H_0$.

**HS interaction matrices.** Recall from the main notes that the interaction is written as $\hat H_I = \sum_i \lambda_i \hat A_i^2$, where each $\hat A_i = \sum_{ab} [A_i]_{ab} c_a^\dagger c_b + a_i$ is a Hermitian fermion bilinear plus a c-number constant $a_i$. After HS decoupling, each vertex contributes a factor $e^{\eta_l \hat V_i}$ where

$$
\hat V_i = \sqrt{-\lambda_i\Delta\tau} \hat A_i
        = \sum_{ab} [V_i]_{ab} c_a^\dagger c_b + v_i,
\quad
V_i = \sqrt{-\lambda_i\Delta\tau} A_i,
\quad
v_i = \sqrt{-\lambda_i\Delta\tau} a_i
$$

For the on-site Hubbard model $V_i$ is a scalar; we precompute the local matrix propagators $e^{\eta_l V_i}$ for every $l \in \lbrace \pm1, \pm2 \rbrace$. The constant $v_i$ is handled through the acceptance ratio (see below).

**HS decoupling channels.** The Hubbard interaction admits two standard decouplings:

**Spin channel.**

$$
U n_{i\uparrow} n_{i\downarrow}
= -\frac{U}{2}(n_{i\uparrow} - n_{i\downarrow})^2
  + \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow})
$$

so after absorbing the chemical potential shift into the single-body Hamiltonian,

$$
\lambda_i = -\frac{U}{2},
\quad
\hat A_i = n_{i\uparrow} - n_{i\downarrow}
=
\begin{bmatrix} c_{i\uparrow}^\dagger & c_{i\downarrow}^\dagger \end{bmatrix}
\begin{bmatrix}1 & \\ & -1\end{bmatrix}
\begin{bmatrix} c_{i\uparrow} \\ c_{i\downarrow} \end{bmatrix}
$$

The single-particle coupling is

$$
[V_i]_{\uparrow\uparrow} = -[V_i]_{\downarrow\downarrow} = \sqrt{-\lambda_i\Delta\tau} = \sqrt{\Delta\tau U/2},
\qquad
v_i = 0.
$$

We store both spin blocks explicitly ($N_\text{f} = 2$).

**$SU(2)$ (charge) channel.**

$$
U n_{i\uparrow} n_{i\downarrow}
= \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow} - 1)^2
  - \frac{U}{2}(n_{i\uparrow} + n_{i\downarrow}) + \frac{U}{2}
$$

so after absorbing the chemical potential shift,

$$
\lambda_i = \frac{U}{2},
\quad
\hat A_i = n_{i\uparrow} + n_{i\downarrow} - 1
=
\begin{bmatrix} c_{i\uparrow}^\dagger & c_{i\downarrow}^\dagger \end{bmatrix}
\begin{bmatrix}1 & \\ & 1\end{bmatrix}
\begin{bmatrix} c_{i\uparrow} \\ c_{i\downarrow} \end{bmatrix}
- 1
$$

The single-particle coupling per stored block is

$$
[V_i]_{\uparrow\uparrow} = [V_i]_{\downarrow\downarrow} = \sqrt{-\lambda_i\Delta\tau} = \sqrt{-\Delta\tau U/2},
\qquad
v_i = -\sqrt{-\lambda_i\Delta\tau}/2.
$$

We store only one fermion block and square the determinant. This HS decoupling is $SU(2)$ symmetric: the decoupled single-particle interaction is invariant under

$$
c_{i\sigma}^\dagger \to c_{i\sigma'}^\dagger [U_i]_{\sigma'\sigma}
$$

for any $SU(2)$ matrix $U_i$.

### Projector parameters and auxiliary field configuration

We want to evaluate the projector partition function

$$
Z_\Theta = \braket{\Psi_L | e^{-(2\Theta + \beta)\hat H} | \Psi_R},
\quad
N_\tau \Delta\tau = 2\Theta + \beta
$$

where $2\Theta+\beta$ is the total projection length and $\beta$ is a window that enables time-displaced measurements. A Monte Carlo configuration is the space–time auxiliary-field array

$$
\mathcal{C} = \{l_{t,i}\}, \quad
t = 1, \dots, N_\tau,  i = 1, \dots, N_I
$$

where $N_I$ is the number of interaction vertices, $l_{t,i} \in \lbrace \pm1, \pm2 \rbrace$ labels the discrete HS field at imaginary-time slice $t$ and interaction vertex $i$ (for the Hubbard model there is one vertex per site, so $N_I = N_s$).

**Trial wavefunction.** We use a right trial state $\ket{\Psi_R}$ and a left trial state $\bra{\Psi_L}$, both Slater determinants built from the same $N \times N_\text{part}$ orbital matrix $P$:

$$
\ket{\Psi_R} = \prod_{j=1}^{N_\text{part}}
\left(\sum_{i=1}^N c_i^\dagger P_{ij}\right) \ket{0}, \quad
P_R = P  (N \times N_\text{part})
$$

$$
\bra{\Psi_L} = \bra{0} \prod_{j=1}^{N_\text{part}}
\left(\sum_{i=1}^N P^\dagger_{ji} c_i\right),
\quad
P_L = P^\dagger  (N_\text{part} \times N)
$$

In practice the orbitals in $P$ are the lowest $N_\text{part}$ eigenvectors of a trial Hamiltonian $T_\text{trial} = T + \delta T_\text{stagger}$, where the small staggered shift lifts degeneracies at the Fermi level (one can refer to the code for an example).

### Discretized propagators and the equal-time Green's function

For a fixed configuration $\mathcal{C}$, the many-body propagator on slice $t$ is the operator

$$
U_t = e^{-\Delta\tau \hat H_0} \prod_i e^{\eta_{l_{t,i}} \hat V_i},
\quad
\hat V_i = \sum_{ab} [V_i]_{ab} c_a^\dagger c_b + v_i
$$

Each $U_t$ is the exponential of a fermion bilinear; it is fully characterised by its $N\times N$ single-particle matrix

$$
B_t = e^{-\Delta\tau T} \prod_i e^{\eta_{l_{t,i}} V_i}
$$

The full many-body propagator is $U_\mathcal{C} = U_{N_\tau} \cdots U_1$, and its matrix representation is $B_\mathcal{C} = B_{N_\tau} \cdots B_1$. The fermion weight of $\mathcal{C}$ is the determinant

$$
W[\mathcal{C}] = \left(\prod_{t,i} \gamma_{l_{t,i}} e^{\eta_{l_{t,i}} v_i}\right)\det\!\left[P_L B_\mathcal{C} P_R\right]
$$

Rather than multiply out $B_\mathcal{C}$ explicitly, we work with the left and right **projected wavefunctions** at an equal-time cut $t$:

$$
\begin{aligned}
B_\rangle(t) &= B_t \cdots B_1 P_R
\quad (N \times N_\text{part}), \\[2pt]
B_\langle(t) &= P_L B_{N_\tau} \cdots B_{t+1}
\quad (N_\text{part} \times N).
\end{aligned}
$$

The equal-time Green's function is defined by the projector expectation value

$$
G_{ij}(t)
\equiv
\frac{
\braket{\Psi_L | U_{N_\tau} \cdots U_{t+1}c_i c_j^\dagger U_t \cdots U_1 | \Psi_R}
}{
\braket{\Psi_L | U_{N_\tau} \cdots U_1 | \Psi_R}
}
$$

Using $B_\rangle(t)$ and $B_\langle(t)$, a standard determinant identity reduces this to the compact matrix expression

$$
G(t) = I - B_\rangle(t)\left[B_\langle(t) B_\rangle(t)\right]^{-1}B_\langle(t)
$$

This formula is the computational backbone of the entire algorithm: local acceptance rate, equal-time observables are all explicitly depends on $G(t)$.

### Local updates and the fast acceptance ratio

For the derivation of this section, please refer to [Assaad and Evertz, Lect. Notes Phys. 739, 277 (2008)](https://link.springer.com/chapter/10.1007/978-3-540-74686-7_10).

If we proposed a local field change and recomputed $W[\mathcal{C}]$ from scratch, each move would cost $O(N_\tau N^3)$ — impossible for all but the smallest systems. The key insight is that we can evaluate the acceptance ratio from the equal-time Green's function at the update location.

When the HS field at $(t,i)$ changes from $l$ to $l'$, the new configuration $\mathcal{C}'$ differs from $\mathcal{C}$ only at that one field. The Metropolis acceptance ratio is

$$
r \equiv \frac{W[\mathcal{C}']}{W[\mathcal{C}]}
   = \frac{\gamma_{l'}}{\gamma_l}
     \cdot
     \frac{e^{\eta_{l'} v_i}}{e^{\eta_l v_i}}
     \cdot
     \frac{\det[P_L B_{\mathcal{C}'} P_R]}{\det[P_L B_{\mathcal{C}} P_R]}
$$

The first factor $\gamma_{l'}/\gamma_l$ is the quadrature weight ratio; the second factor comes from the constant shift $v_i$; the third is the determinant ratio. To evaluate the determinant ratio efficiently we define the local update matrix

$$
\Delta = e^{\eta_{l'} V_i} e^{-\eta_l V_i} - I
$$

which is nonzero only on a few orbitals touched by that vertex. The determinant ratio then reduces to

$$
\frac{\det[P_L B_{\mathcal{C}'} P_R]}{\det[P_L B_{\mathcal{C}} P_R]}
= \det\!\left[I + (I - G)\Delta\right]
$$

where $G$ is the equal-time Green's function evaluated just before the updated vertex. For the on-site Hubbard model the vertex acts only on site $i$, so $\Delta$ is nonzero only at the spin-diagonal entries $\Delta_{i\sigma,i\sigma}$ ($\sigma = \uparrow,\downarrow$). The determinant factorises over spin flavours, giving per flavour

$$
1 + \Delta_{i\sigma,i\sigma}\left(1 - G_{i\sigma,i\sigma}\right)
$$

The constant-shift factor is $\exp[v_i(\eta_{l'} - \eta_l)]$ per stored fermion flavor. Taking the product over all flavors, the full local ratio is

$$
r = \frac{\gamma_{l'}}{\gamma_l}
    \prod_{\sigma=1}^{N_\text{f}}
    e^{v_i(\eta_{l'} - \eta_l)}
    \left[ 1 + \Delta_{i\sigma,i\sigma}\left(1 - G_{i\sigma,i\sigma}\right) \right]
$$

For the $SU(2)$ channel a single flavor block is stored but represents both spin species; the determinant is squared accordingly.

**Sherman–Morrison update.** When a proposal is accepted we update $G$ via the rank-1 Sherman–Morrison formula. Let $k$ be the orbital touched by the vertex. The update reads

$$
G_{ij} \leftarrow G_{ij} + \frac{\Delta_{kk}}{1 + \Delta_{kk}(1 - G_{kk})}  G_{ik}  (\delta_{kj} - G_{kj})
$$

which costs $O(N^2)$ rather than the $O(N^3)$ of recomputing $G$ from scratch. Together with the $O(1)$ ratio evaluation, a full sweep costs $O(N_\tau N_s N^2)$ — the per-slice Green's-function propagation and per-site update are both $O(N^2)$, and the expensive $O(N^3)$ determinant inversions are confined to stabilization boundaries.

### Sweep structure

We alternate forward ($t = 1 \to N_\tau$) and backward ($t = N_\tau \to 1$) sweeps. In the forward direction we first advance the Green's function through the kinetic factor,

$$
G \leftarrow e^{-\Delta\tau T} G e^{+\Delta\tau T}
$$

and then through each local interaction vertex on that slice,

$$
G \leftarrow e^{\eta_l V_i} G e^{-\eta_l V_i}
$$

After propagating through a vertex we propose a new HS field value, compute the ratio $r$, and accept with probability $\min(1, \lvert \mathrm{Re}\,r \rvert)$. If accepted we update the field and correct $G$ via Sherman–Morrison. We measure equal-time observables at the symmetric midpoint $t = N_\tau/2$, and rebuild $G$ from cached projectors whenever $t$ is a multiple of $n_\text{stab}$.

The backward sweep traverses the slices in reverse. The vertex loop runs in reverse order, and the kinetic propagation uses the inverse factor $G \leftarrow e^{+\Delta\tau T} G e^{-\Delta\tau T}$. We collect time-displaced observables only on the backward sweep.

### Measurements

**Equal-time observables.** At $t = N_\tau/2$ we measure observables using the equal-time Green's function $G$. For the HS-decoupled system the auxiliary fields are fixed, so fermion operators obey Wick's theorem: any multi-particle expectation value factorizes into a sum over products of single-particle contractions. Assuming $\braket{c^\dagger c^\dagger} = \braket{cc} = 0$ (no pairing), the only nonzero contractions are $\braket{c_i^\dagger c_j} = \delta_{ij} - G_{ji}$ and $\braket{c_i c_j^\dagger} = G_{ij}$.

For example, a four-fermion expectation (under a single configuration $\mathcal{C}$) equals

$$
\braket{c_i^\dagger c_j c_k^\dagger c_l}
= \braket{c_i^\dagger c_j} \braket{c_k^\dagger c_l}
  + \braket{c_i^\dagger c_l} \braket{c_j c_k^\dagger}
$$

Two important observables are:

- **Kinetic energy:**
  $$
  \braket{\hat H_0} = \sum_{ij} T_{ij} \braket{c_i^\dagger c_j} = \text{Tr}[T (I - G)]
  $$
- **Potential energy:** $\braket{\hat H_I} = U \sum_i \braket{n_{i\uparrow} n_{i\downarrow}}$. Applying Wick's theorem,
  $$
  \braket{n_{i\uparrow} n_{i\downarrow}}
  = \braket{c_{i\uparrow}^\dagger c_{i\uparrow}} \braket{c_{i\downarrow}^\dagger c_{i\downarrow}}
  = (1 - G_{i\uparrow,i\uparrow})(1 - G_{i\downarrow,i\downarrow})
  $$

All observables are accumulated weighted by the current sign. The measurement corresponds to the symmetric estimator

$$
\braket{\hat O}
=
\frac{
\braket{\Psi_L | e^{-(\Theta + \beta/2)\hat H}\hat Oe^{-(\Theta + \beta/2)\hat H} | \Psi_R}
}{
\braket{\Psi_L | e^{-(2\Theta + \beta)\hat H} | \Psi_R}
}
$$

which converges to the ground-state expectation as $\Theta \to \infty$.

**Time-displaced Green's functions**

$$
G_{ij}(\tau) \equiv \braket{c_i(\tau) c_j^\dagger(0)},
\quad
G_{ij}(-\tau) \equiv -\braket{c_j^\dagger(\tau) c_i(0)}
$$

for $\tau = 0, \Delta\tau, \dots, \beta$. The measurement is anchored at $t_0 = \Theta/\Delta\tau$. At each stabilization boundary during the propagation we reconstruct the equal-time $G$ from cached projectors to prevent roundoff error from accumulating over long imaginary-time separations.

### Sign, reweighting, and binning

The fermion weight $W[\mathcal{C}]$ is not guaranteed to be real and positive. As discussed in the main notes, we sample with respect to $\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert$ and recover physical expectation values by reweighting with

$$
\text{sign}[\mathcal{C}] = \frac{W[\mathcal{C}]}{\lvert \mathrm{Re}\,W[\mathcal{C}] \rvert}
$$

We update the sign incrementally: for each accepted move, $\text{sign} \leftarrow \text{sign} \times \mathrm{Re}\,r / \lvert \mathrm{Re}\,r \rvert$. The Metropolis acceptance probability is $A = \min(1, \lvert \mathrm{Re}\,r \rvert)$.

We divide the production run into bins of $N_\text{sweep}$ forward–backward sweeps each. Within a bin, signed observables are accumulated and normalised by the number of measurements. Independent MPI ranks run separate Markov chains with different random seeds; we reduce bin averages across ranks. Final error analysis uses jackknife resampling — plain jackknife for the average sign, and jackknife ratios for sign-reweighted observables:

$$
\braket{\hat O}
=
\frac{
\braket{\text{sign}\cdot \hat O}_\text{MC}
}{
\braket{\text{sign}}_\text{MC}
}
$$

### Numerical stabilization

Products of many non-unitary matrices are exponentially ill-conditioned: the smallest singular values of $B_t \cdots B_1$ decay as $e^{-\text{const}\times t}$, quickly destroying numerical precision. The standard remedy is periodic re-orthonormalization by a thin QR decomposition, keeping only the $Q$ factor. We denote the number of slices between stabilizations by $n_\text{stab}$.

The stabilization operates at three levels. At startup we build the full forward and backward projector chains and cache the $Q$ factors at every $n_\text{stab}$-th slice. At the beginning of each bin, and whenever we cross a stabilization boundary during a sweep, we reconstruct $G$ exactly from the nearest cached stabilization point — this costs $O(N_\tau N^3 / n_\text{stab})$ per bin. At each stabilization boundary during a sweep we refresh the cached projector to reflect accumulated Monte Carlo updates, recompute $G$, and abort if the discrepancy with the propagated $G$ exceeds $10^{-3}$.

## Minimal implementation
A minimal DQMC program for the 1D/2D Hubbard model written in Julia is under the [hubbard-afqmc](https://github.com/haoran0115/haoran0115.github.io/tree/main/notes/dqmc_tut/hubbard-afqmc) directory
* `hubbard-afqmc.jl`: utility functions
* `main.jl`: main QMC loop
* `parameters.jl`: simulation parameters
* `seeds.jl`: random seed
* `analysis.jl`: data analysis

How to run the code: modify the simulation parameters in `parameters.jl` and run with
```bash
# with single Markov chain
julia main.jl

# with N Markov chains, e.g. N = 4 for the following command
mpiexecjl -np 4 julia main.jl
```

Simulation results will be written under `data/` and analyzed results will be written under `analysis/`.

## Reference and useful resources

Tutorials
* [TowardQMC](https://www.youtube.com/playlist?list=PLheYERt_Ks3wHU_Fa7MVGH5I31pJf4SLH): video introduction for finite-$T$ formalism
* [DQMC](https://quantummc.xyz/teaching/dqmc/): concise introduction for finite-$T$ DQMC using Hubbard model as an example, written in Chinese

Technical walk-through
* [World line and determinantal Quantum Monte
Carlo methods for spins, phonons, and
electrons](https://pawn.physik.uni-wuerzburg.de/~assaad/Reprints/assaad_evertz.pdf): these lecture notes form the basis of the DQMC part of this tutorial, and the minimal implementation closely follows the second half (DQMC) of these notes

Software packages
* [ALF](https://git.physik.uni-wuerzburg.de/ALF/ALF): very comprehensive DQMC/AFQMC package in Fortran, with complete documentation, minimal tutorials, and [workshop recordings&tutorials](https://git.physik.uni-wuerzburg.de/ALF/ALF_Tutorial/-/tree/master/Presentations?ref_type=heads)
* [SmoQyDQMC](https://smoqysuite.github.io/SmoQyDQMC.jl/stable/): HMC (a determinant-free QMC) for electron-electron and electron-phonon problems
