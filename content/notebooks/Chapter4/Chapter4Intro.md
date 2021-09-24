# Chapter 4 Controllability , observability and stability

In the last chapter we obtained a formula for the solution of a linear system of first order differential equations of the form 

$$
 \frac{\md\mathbf{x}}{\md t} = A\mathbf{x} + \mathbf{f}(t)\quad : \quad \mathbf{x}(0) = \mathbf{x}_0,
$$ (consys)

*i.e.*,

$$
\mathbf{x}(t) = e^{At}\mathbf{x}_0+\int_0^te^{A(t-\tau)}\mathbf{f}(\tau)\mathrm{d}\tau.
$$

Suppose that we now rewrite the *forcing effect* in a form that shows explicitly how many *control functions* can be used to try to influence the behaviour of the system,

$$
 \text{i.e., } \mathbf{f}(t) = B\mathbf{u}(t),
$$

where $B$ is a constant $n\times p$ matrix and $\mathbf{u}(t)$ is a column $p$-vector consisting of these `control functions'. In these circumstances the solution of a system of the above form at time $t=t_f$ is given by

$$
\mathbf{x}(t_f) &= e^{At_f}\mathbf{x}_0+\int_0^{t_f}e^{A(t_f-\tau)}B\mathbf{u}(\tau)\md\tau,\\ \label{eqn:controlode}
&= e^{At_f}\mathbf{x}_0+\int_0^{t_f}e^{A\tau}B\mathbf{u}(t_f-\tau)\md\tau.
$$

In this chapter we shall derive conditions under which systems of the form in {eq}`consys` can be controlled.
