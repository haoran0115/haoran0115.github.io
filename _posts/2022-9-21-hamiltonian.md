---
title: 'An interesting property of Hamitonian system'
tags:
  - Classical mechanics
  - Differential equations
---
An interesting problem appeared in a book about classical mechanics.[^fn]
> Problem. Show that in an hamiltonian system it is impossible to have asymptotically stable equilibrium positions and asymptotically stable limit cycles in the phase space.

*Ideas:* first consider a system of odes
$$
\dot{\mathbf{x}} = f(\mathbf{x})
$$
In a neighborhood around an asymptotically stable equilibrium position $\mathbf{x}_0$, the equation could be
written in expansion
$$
\dot{\mathbf{x}} = A\mathbf{x} + \mathbf{o}(\|\mathbf{x}\|^2)
$$
where the matrix $A$ with all eigenvalue negative, and $\mathbf{o}$ is the remainder.
To prove a hamiltonian system $H(\mathbf{q}, \dot{\mathbf{q}}, t)$ cannot have any asymptotically stable limit point, we just show that the matrix $A$ for the equation
$$
\dot{\mathbf{q}} = A\mathbf{q} + \mathbf{o}(\|\mathbf{q}\|^2)
$$

The case for the limit cycle is similar.



[^fn]: Arnol'd, V. I. (2013). Mathematical methods of classical mechanics (Vol. 60). Springer Science & Business Media.

