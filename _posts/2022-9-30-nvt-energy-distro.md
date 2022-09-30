---
title: Energy distribution under NVT ensemble
tag: 
    - Statistical mechanics
    - Physics
---


Under the NVT ensemble, the probability density of a particular state is given by the famous Boltzmann distribution. Note that $$Q(N,V,T)$$ is the partition function, $\mathbf{x}$ denotes the phase space coordinate, $\beta=1/kT$, and $\mathcal{H}(\mathbf{x})$ is the Hamiltonian of the system.

$$
f(\mathbf{x}) = \frac{e^{-\beta\mathcal{H}(\mathbf{x})}}{Q(N,V,T)},\
Q(N,V,T) = \int \mathrm{d}\mathbf{x}\ e^{-\beta\mathcal{H}(\mathbf{x})}
$$

In order to derive the probability density $P(E)$ under NVT ensemble, one way is to perform an integration over $\mathbf{x}$ since we have the relation $E(\mathbf{x}) = \mathcal{H}(\mathbf{x})$. Some people may add a normalization constant $C_N=1/(h^{3N}N!)$ to $f(\mathbf{x})$, however, I dropped this constant for convenience. Please also note that $\Omega(N,V,E)$ is the partition function of the NVE (microcanonical) ensemble.

$$
\begin{align}
\begin{split}
    P(E) &= \int_{\mathcal{H}(x)=E}\mathrm{d}\mathbf{x}\ f(\mathbf{x})\\
    &= \int\mathrm{d}\mathbf{x}\ f(\mathbf{x})\delta(\mathcal{H}(\mathbf{x})-E)\\
    &= \int\mathrm{d}\mathbf{x}\frac{e^{-\beta\mathcal{H}(\mathbf{x})}}{Q(N,V,T)}\delta(\mathcal{H}(\mathbf{x})-E)\\
    &= \int\mathrm{d}\mathbf{x}\ \frac{e^{-\beta E}}{Q(N,V,T)}\delta(\mathcal{H}(\mathbf{x})-E)\\
    &= \frac{e^{-\beta E}}{Q(N,V,T)}\int\mathrm{d}\mathbf{x}\ \delta(\mathcal{H}(\mathbf{x})-E)\\
    &= \frac{\Omega(N,V,E)}{Q(N,V,T)}e^{-\beta E}
\end{split}
\end{align}
$$

One application of this conclusion is provided as the following problem.
> Problem. Please describe the relationship between the probability and energy of state in NVT ensemble for a system contains $N$ one-dimensional harmonic oscillators.


*Solution.* Since the partition function for $N$ one-dimensional harmonic oscillators under the NVE ensemble is in the form of $\Omega(N, E)=f(N)E^{N-1}$. Therefore, under the NVT ensemble, we would have

$$
P(E)\propto \Omega(N,E)e^{-\beta E} \propto E^{N-1}e^{-\beta E}
$$


