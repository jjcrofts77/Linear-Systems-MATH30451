# Chapter 5 System and transfer functions

Consider a system whose *input-state* and *state-output* relationship are given by

$$
 \frac{\mathrm{d}\mathbf{x}}{\mathrm{d} t} &= A\mathbf{x}+B\mathbf{u},\\
 \mathbf{y}&=C\mathbf{x},
$$

where $A, B$ and $C$ are $n\times n, n\times p$ and $l\times n$ matrices, respectively.

Then, assuming that $\ds\mathbf{x}(0)=0$, taking Laplace transforms gives

$$
 s\bar{\mathbf{x}} &= A\bar{\mathbf{x}}+B\bar{\mathbf{u}},\\
 \bar{\mathbf{y}} &=C\bar{\mathbf{x}}.
$$

Rearranging the first of these two equations to give

$$
 (sI_n-A)\bar{\mathbf{x}} = B\bar{\mathbf{u}},
$$

we see that the situation can be represented in matrix form as

$$
 \begin{pmatrix} sI_n-A&B\\-C&0\end{pmatrix}\begin{pmatrix}\bar{\mathbf{x}}\\-\bar{\mathbf{u}}\end{pmatrix} = \begin{pmatrix}0\\-\bar{\mathbf{y}}\end{pmatrix}.
$$

Here, the $(n+l)\times(n+p)$ matrix 

$$
 P(s) = \begin{pmatrix} sI_n-A&B\\-C&0\end{pmatrix},
$$

is called the *state-space system matrix* of the given system.

Now, from the rearranged equation involving $\bar{\mathbf{x}}$ and $\bar{\mathbf{u}}$,
we see that

$$
 \bar{\mathbf{x}} = (sI_n-A)^{-1}B\bar{\mathbf{u}},
$$

and hence that

$$
 \bar{\mathbf{y}} = C(sI_n-A)^{-1}B\bar{\mathbf{u}} = G(s)\bar{\mathbf{u}}.
$$

Here we have used the $l\times p$ matrix 

$$
 G(s) = C(sI_n-A)^{-1}B,
$$

which is called the *transfer function matrix* of the system. Knowledge
of this matrix allows us to convert the Laplace transform of the input, 
$\bar{\mathbf{u}}$, into that of the output, $\bar{\mathbf{y}}$.

<br>

**Example 5.1** The system

$$
 \frac{\md}{\md t}\bmc x_1\\x_2\\x_3\end{pmatrix} &= \begin{pmatrix}6&-2&-5\\-1&-1&1\\8&-2&-7\end{pmatrix}\begin{pmatrix} x_1\\x_2\\x_3\end{pmatrix}+\begin{pmatrix}2\\0\\3\end{pmatrix} u,\\
 y &= [1,-2,0]\begin{pmatrix} x_1\\x_2\\x_3\end{pmatrix},
$$

gives rise to the state-space matrix

$$
 P(s) = \begin{pmatrix} s-6&2&5&|&2\\1&s+1&-1&|&0\\-8&2&s+7&|&3\\-&-&-&|&-\\-1&2&0&|&0\end{pmatrix},
$$

and the transfer function (matrix)

$$
 G(s) &= \begin{pmatrix}1&-2&0\end{pmatrix}\begin{pmatrix} s-6&2&5\\1&s+1&-1\\-8&2&s+7\end{pmatrix}^{-1}\begin{pmatrix}2\\0\\3\end{pmatrix}\\
      &=\begin{pmatrix} 1&-2&0\end{pmatrix}\frac{1}{\Delta(s)}\begin{pmatrix} s^2+8s+9&*&-5s-7\\1-s&*&s-1\\ *&*&*\end{pmatrix}
      \begin{pmatrix}2\\0\\3\end{pmatrix},
$$

where

$$
\left|\begin{matrix}s-6&2&5\\1&s+1&-1\\-8&2&s+7\end{matrix}\right|
&=(s-6)(s^2+8s+9)+2(1-s)+5(8s+10)\\
&=s^3+2s^2-s-2=(s+2)(s^2-1)\\
&=(s+2)(s+1)(s-1).
$$

Then 

$$
 G(s) &= \frac{1}{\Delta(s)}\begin{pmatrix}1&-2&0\end{pmatrix}\begin{pmatrix}2s^2+s-3\\s-1\\ *\end{pmatrix} = \frac{2s^2-s-1}{(s+2)(s+1)(s-1)},\\
      &=\frac{2s+1}{(s+2)(s+1)}.
$$

Now, suppose that the system has input $\displaystyle u=2e^{-3t}$, for which $\displaystyle\bar{u}=\frac{2}{s+3}$. Then

$$
 \bar{y} &= \frac{2s+1}{(s+1)(s+2)}\cdot\frac{2}{s+3} = \frac{4s+2}{(s+1)(s+2)(s+3)},\\
         &= -\frac{1}{s+1}+\frac{6}{s+2}-\frac{5}{s+3}.
$$

Hence $\displaystyle y=-e^{-t}+6e^{-2t}-5e^{-2t}$.

Now, while it is clear that a given system or state-space system matrix can
give rise to only one corresponding transfer function matrix, the reverse is
not true. A given transfer function matrix can be the result of several 
state-space matrices or systems, *i.e.* it has many possible realisations.

<br>

**Example 5.2** The above transfer function can be realised by any 
of the following state-space system matrices

$$
 \text{(i) } \begin{pmatrix} s-1&0&0&|&0\\0&s+1&0&|&1\\3&0&s+2&|&1\\-&-&-&|&-\\-4&1&-3&|&0\end{pmatrix}
 &\qquad \text{(ii) } \begin{pmatrix} s&0&-1&|&0\\-3&s&1&|&1\\2&0&s+3&|&1\\-&-&-&|&-\\-1&0&-2&|&0\end{pmatrix}\\
 \text{(iii) } \begin{pmatrix} s+1&0&|&1\\0&s+2&|&1\\-&-&|&-\\1&-3&|&0\end{pmatrix} &\qquad \text{(iv) }
 \begin{pmatrix} s&-1&|&0\\2&s+3&|&1\\-&-&|&-\\-1&-2&|&0\end{pmatrix}.
$$

Note that the system matrix in (iii) is the result of removing all reference to the
*input-decoupling zero* $s=1$ from that in (i), while that in (iv) is the 
result of similarly removing reference to the *output-decoupling zero* 
$s=0$ from that in (ii). The system matrices in (iii) and (iv) may be thought of as giving more efficient realisations of the corresponding transfer function matrix since they are smaller than those in (i) and (ii). They do, in fact, provide *minimal realisations* of the transfer function involved.

System and transfer function matrices can also arise from systems of difference 
equations. In this case the Laplace transform parameter $s$ is replaced by the
z-transform parameter $z$.

<br>

**Example 5.3** Taking z-transforms of the system of difference equations

$$
 \begin{pmatrix} x^{(1)}_{k+1}\\x^{(2)}_{k+1}\\x^{(3)}_{k+1}\end{pmatrix} = 
 \begin{pmatrix}6&-2&-5\\-1&-1&1\\8&-2&-7\end{pmatrix}
 \begin{pmatrix} x^{(1)}_k\\x^{(2)}_k\\x^{(3)}_k\end{pmatrix}+
 \begin{pmatrix}2\\0\\3\end{pmatrix} u_k,
$$

gives rise to the system matrix

$$
 P(z) = \begin{pmatrix} z-6&2&5&|&2\\1&z+1&-1&|&0\\-8&2&z+7&|&3\\-&-&-&|&-\\-1&2&0&|&0\end{pmatrix},
$$

and the transfer function 

$$
 G(z) = \frac{2z+1}{(z+1)(z+2)}.
$$

Then, if the above system has input $\displaystyle u_k = 2\left(-3\right)^k$, for which 
$\displaystyle \bar{u}=\frac{2z}{z+3}$, the response of the system has z-transform

$$
 \bar{y} &= \frac{2z+1}{(z+1)(z+2)}\times\frac{2z}{z+3} = \frac{(4z+2)z}{(z+1)(z+2)(z+3)},\\
         &=\left(\frac{-1}{z+1}+\frac{6}{z+2}-\frac{5}{z+3}\right)z,
$$

*i.e.* $\displaystyle y_k = (-1)^k+6(-2)^k-5(-3)^k$.
