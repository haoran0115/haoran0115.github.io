---
title: "Liouville&#39;s formula"
tags:
  - Differential equations
---
> Problem. Prove Liouville&#39;s formula $$W=W_0 e^{\int \mathrm{tr}~ A dt}$$ for the Wronskian determinant of the linear system $$\dot{\mathbf{x}} = A(t)\mathbf{x}$$.

*Proof.* Let $$\mathbf{x}(t)\in\mathbb{R}^n$$ and $$A(t)\in\mathbb{R}^{n\times n}$$. The solution to the system $$\dot{\mathbf{x}} = A(t)\mathbf{x}$$ is in the following form.

$$
\mathbf{x}(t) = e^{\int_0^t A(s)~ ds}\mathbf{x}(0)
$$

Since $$e^{\int_0^t A(s)~ ds}$$ is invertible, there would have $n$ linearly independent solutions $\mathbf{x}_1,\dots,\mathbf{x}_n$. Hence, we can write the Wronskian in the following form.

$$
W(t) = \begin{bmatrix}
    \mathbf{x}_1(t) & \cdots & \mathbf{x}_n(t)
\end{bmatrix} = 
e^{\int_0^t A(s)~ ds}
\begin{bmatrix}
    \mathbf{x}_1(0) & \cdots & \mathbf{x}_n(0)
\end{bmatrix} = 
e^{\int_0^t A(s)~ ds} W(0)
$$

Therefore $$\det W(t) = \det e^{\int_0^t A(s)~ ds}\det W(0)$$.

*Claim.* $$\det e^A = e^{\mathrm{tr}~ A}$$

> *Proof of Claim.* Decompose $$A=UTU^{-1}$$, where $$T$$ is an upper triangular matrix. Then we can > show that $$\det e^{B}=e^{T}$$ by the following derivation.
> 
> $$
\begin{aligned}
e^{A} &= \sum_n \frac{1}{n!} A^n
= \sum_n\frac{1}{n!} (UTU^{-1})^n
= \sum_n\frac{1}{n!} UT^n U^{-1}
= U \left(\sum_n\frac{1}{n!}\right) U^{-1}
= U e^T U^{-1}\\
\Rightarrow \det e^A &= \det U \det e^T \det U^{-1} = \det e^T    
\end{aligned}$$
> 
> Since $T$ is upper-triangular, then $$e^T$$ is also upper-triangular. Due to the property of upper-triangular matrices where $$[e^T]_{ii} = e^{T_{ii}}$$, and $$\det T = \prod_i T_ {ii}$$, we can prove that $$\det e^T = e^{\mathrm{tr}~ T}$$.
> 
> $$
\det e^T = \prod_i [e^T]_{ii} = \prod_i e^{T_{ii}} = e^{\sum_i T_{ii}} = e^{\mathrm{tr}~ T}$$
> 
> Since $A$ and $T$ has the same set of eigenvalues, we have $\mathrm{tr}~ A=\mathrm{tr}~ T$ Therefore 
> 
> $$
\det e^A = \det e^T = e^{\mathrm{tr}~ T} = e^{\mathrm{tr}~ A}$$

By the claim we just proved, one can show that $$\det W(t) = e^{\int_0^t \mathrm{tr}~ A(s)~ ds}\det W(0)$$.

$$
\det W(t) = \det W(t) =  e^{\mathrm{tr}\int_0^t A(s)~ ds} \cdot \det W(0)
=  e^{\int_0^t \mathrm{tr}~ A(s)~ ds} \cdot \det W(0)
$$





