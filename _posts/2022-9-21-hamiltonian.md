---
title: 'An interesting property of Hamitonian system'
tags:
  - Classical mechanics
  - Differential equations
---
An interesting problem appeared in a book about classical mechanics.[^fn]
> Problem. Show that in an hamiltonian system it is impossible to have asymptotically stable equilibrium positions and asymptotically stable limit cycles in the phase space.

*Ideas:* first consider a system of ODEs
\\[
\dot{\mathbf{x}} = f(\mathbf{x})
\\]
In a neighborhood around an asymptotically stable equilibrium position $\mathbf{x}_0$, the equation could be
written in expansion
\\[
\dot{\mathbf{x}} = A(\mathbf{x}-\mathbf{x}_0) + O(\|\mathbf{x}-\mathbf{x}_0\|^2)
\\]
where the matrix $A$ with all eigenvalue negative, and $O$ is the remainder (which is much small when $\|\mathbf{x}-\mathbf{x}_0\|$ is small).
Then, one way to prove a hamiltonian system $H(\mathbf{q}, \dot{\mathbf{q}}, t)$ cannot have any asymptotically stable limit point is to show that the matrix $A$ for the equation
\\[
  \dot{\mathbf{q}} = A(\mathbf{q}-\mathbf{q}_0) + O(\|\mathbf{q}-\mathbf{q}_0\|^2)
\\]
cannot have all its eigenvalues negative/have negative real part.

The case for the limit cycle would be much more complicated, but the idea is similar.

*Proof.* Let $\mathbf{x}=[q_1, \dots, q_n, p_1, \dots, p_n]^T$ denotes the phase space coordinate. Then, a Hamiltonian system satisfies
\\[
  \dot{q}_i = \frac{\partial H}{\partial p_i},\
  \dot{p}_i = -\frac{\partial H}{\partial p_i}
\\]
Then the equation could be written in the form of
\\[
  \dot{\mathbf{x}} = 
  \begin{bmatrix}
  \dot{q}_1\\\\ \cdots\\ \dot{q}_n \\ \dot{p}_1\\ \cdots \\ \dot{p}_n
  \end{bmatrix}
  = \begin{bmatrix}
    \frac{\partial H}{\partial p_1}\\ \cdots\\
    \frac{\partial H}{\partial p_n}\\ 
    -\frac{\partial H}{\partial q_1}\\ \cdots\\
    -\frac{\partial H}{\partial q_n}
  \end{bmatrix}
\\]



[^fn]: Arnol'd, V. I. (2013). Mathematical methods of classical mechanics (Vol. 60). Springer Science & Business Media.

