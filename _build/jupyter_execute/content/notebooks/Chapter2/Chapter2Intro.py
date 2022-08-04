#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/jjcrofts77/Linear-Systems-MATH30451/blob/main/content/notebooks/Chapter2/Chapter2Intro.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# In[ ]:





# # Chapter 2 Similarity Transformations and the Jordan Canonical Form
# 
# Recall that if $A$ is a real symmetric matrix, and if we choose $P$ to be the matrix whose columns are the orthonormal eigenvectors of $A$, then the matrix $D=P^{-1}AP$ is diagonal, with non-zero elements given by the eigenvalues of $A$. 
# 
# ````{margin} Eigen-decomposition
# ```{figure} ../../images/EigDecomp2.png
# ```
# ````
# 
# Importantly, we can use the above decompostion to solve systems of linear equations such as those seen in Lecture 1. For example, using similarity transformations it is straightforward to map a coupled set of ODEs, which are difficult to solve, to an uncoupled set, which can be solved using standard methods, *i.e.*
# 
# $$
#  \frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = A\mathbf{x} \longrightarrow \frac{\mathrm{d}\mathbf{y}}{\mathrm{d}t} = D\mathbf{y},
# $$
# 
# or 
# 
# $$
# \frac{\mathrm{d}}{\mathrm{d}t}\begin{pmatrix} x_1\\x_2\\\vdots\\x_n\end{pmatrix} = \begin{pmatrix} a_{11}&a_{12}&\cdots&a_{1n}\\a_{21}&a_{22}&\cdots &a_{2n}\\\vdots&\vdots&\cdots&\vdots\\a_{n1}&a_{n,2}&\cdots&a_{nn}\end{pmatrix}\begin{pmatrix} x_1\\x_2\\\vdots\\x_n\end{pmatrix}
# \longrightarrow \frac{\mathrm{d}}{\mathrm{d}t}\begin{pmatrix} y_1\\y_2\\\vdots\\y_n\end{pmatrix} =\begin{pmatrix} \lambda_1&&&\\&\lambda_2&&\\&&\ddots&\\&&&\lambda_n \end{pmatrix}\begin{pmatrix} y_1\\y_2\\\vdots\\y_n\end{pmatrix}.
# $$
# 
# For example, the first $x$-equation in the above is given by
# 
# $$
#  \frac{\mathrm{d}x_1}{\mathrm{d}t} = a_{11}x_1+a_{12}x_2
# + \cdots +a_{1n}x_n,
# $$
# 
# whereas the corresponding $y$-equation is
# 
# $$
# \frac{\mathrm{d}y_1}{\mathrm{d}t} = \lambda_1 y_1,
# $$
# 
# which is seperable and so easy to solve.
# 
# <br>
# 
# **Example 2.0.1** Solve the differential equation
# 
# $$
# \mathbf{x}' = A\mathbf{x} = \begin{pmatrix} -3 &1\\1&-3\end{pmatrix}\mathbf{x}
# $$
# 
# via the method of similarity transformations.
# 
# **Solution**
# 
# ```{toggle}
# To start we must determine the eigenvalues and eigenvectors of $A$. The characteristic equation is given by
# 
# $$
# \det(A-\lambda I_2) = \left|\begin{array}{cc}
# -3-\lambda&1\\1&-3-\lambda
# \end{array}\right| = 
# \lambda^2+6\lambda+8 = 0.
# $$
# 
# This gives eigenvalues $\lambda_1=-2$ and $\lambda_2=-4$. Eigenvectors are obtained from 
# 
# $$
# \begin{pmatrix}-1&1\\1&-1\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = 0 \quad \text{or}\quad \begin{pmatrix}-1&~1\\0&0\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = 0 \implies \mathbf{u}_1 = \begin{pmatrix} 1\\1\end{pmatrix} 
# $$
# 
# and
# 
# $$
# \begin{pmatrix}~1&~1\\1&1\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = 0 \quad \text{or}\quad \begin{pmatrix} ~1&~1\\0&0\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = 0 \implies \mathbf{u}_2 = \begin{pmatrix} -1\\1\end{pmatrix}. 
# $$
# 
# Thus
# 
# $$
# P = \begin{pmatrix} \mathbf{u}_1& \mathbf{u}_2\end{pmatrix} = \begin{pmatrix} ~1&-1\\1&1\end{pmatrix} \quad \text{and} \quad P^{-1} = \frac{1}{2}\begin{pmatrix} 1~&1\\-1~&1\end{pmatrix}.
# $$
# 
# By introducing the new variable $\displaystyle\mathbf{y} = P^{-1}\mathbf{x}$ we can rewrite the system of ODEs above as
# 
# $$
#  \frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} &= PDP^{-1}\mathbf{x}\\
# \iff \frac{\mathrm{d}P^{-1}\mathbf{x}}{\mathrm{d}t} &= DP^{-1}\mathbf{x} \qquad&\text{(multiply by $P^{-1}$ on left)}\\
# \iff \frac{\mathrm{d}\mathbf{y}}{\mathrm{d}t} &= D\mathbf{y} \qquad &\text{(using $\mathbf{y} = P^{-1}\mathbf{x}$)}
# $$
# 
# or
# 
# $$
# \frac{\mathrm{d}}{\mathrm{d}t}\begin{pmatrix} y_1\\y_2\end{pmatrix} = \begin{pmatrix} -2&0\\0&-4\end{pmatrix}\begin{pmatrix} y_1\\ y_2\end{pmatrix}.
# $$
# 
# But the uncoupled system above is just a pair of 1st order, seperable ODEs, which we learnt how to solve in first year:
# 
# $$
# \frac{\mathrm{d}y_1}{\mathrm{d}t} =-2y_1 \implies y_1 = Ae^{-2t} \quad\text{ and }\quad \frac{\mathrm{d}y_2}{\mathrm{d}t} =-4y_2 \implies y_2 = Be^{-4t}.
# $$
# 
# Therefore 
# 
# $$
# \mathbf{y} = \begin{pmatrix} Ae^{-2t}\\ Be^{-4t}\end{pmatrix} \implies \mathbf{x} = P\mathbf{y} = \begin{pmatrix} ~1 &-1\\1&1\end{pmatrix}\begin{pmatrix} Ae^{-2t}\\Be^{-4t}\end{pmatrix} = \begin{pmatrix} Ae^{-2t}-Be^{-4t}\\Ae^{-2t}+Be^{-4t}\end{pmatrix}
# $$
# 
# Note that we can rewrite the above solution for $\mathbf{x}(t)$ using the eigenvectors of $A$:
# 
# $$
# \mathbf{x}(t) = A\begin{pmatrix}1\\1\end{pmatrix} e^{-2t} + B\begin{pmatrix} -1\\1\end{pmatrix} e^{-4t}.
# $$
# ```
# <br>
# 
# **Question: What happens when $\mathbf{A}$ is not symmetric?**
# 
# <br>
# Consider, for example,
# 
# $$
# A = \begin{pmatrix} 1 &-1\\1&3\end{pmatrix}.
# $$
# 
# In order to diagonalise $A$ we need to 
# 
# 1. Compute the eigenvalues $\displaystyle D = \begin{pmatrix} \lambda_1&0\\ 0&\lambda_2 \end{pmatrix}$; and
# 
# 2. Compute the eigenvectors $\displaystyle P = \begin{pmatrix} \mathbf{u}_1&\mathbf{u}_2\end{pmatrix}$. 
# 
# The *characterisctic polynomial* of $A$ is given by
# 
# $$
#  \chi_A(t) &= \left|tI_2-A\right|\\
#  &=\left|\begin{array}{cc}
#  t-1&1\\-1&t-3
# \end{array}  \right|\\
# &= (t-1)(t-3)+1\\
# &=t^2-4t+4 = (t-2)^2.
# $$
# 
# Therefore, $\displaystyle\lambda_1=\lambda_2 = 2$, *i.e.* 2 is the only eigenvalue of $A$; we say that the eigenvalue has *algebraic multiplicity* of 2.
# 
# Next we want to compute the *eigenspace* corresponding to $\lambda=2$. To do this we solve
# 
# $$
# \left(\lambda I_2-A\right)\mathbf{u} = \left(2I_2-A\right)\mathbf{u} = 0
# $$
# 
# or
# 
# $$
# \begin{pmatrix} 1&1\\-1&-1\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = \mathbf{0}
# $$
# 
# or
# 
# $$
# \begin{pmatrix} 1&1\\~0&~0\end{pmatrix}\begin{pmatrix} u_1\\u_2\end{pmatrix} = \mathbf{0}.
# $$
# 
# The eigenspace is thus given by
# 
# $$
# \left\{ r\begin{pmatrix} -1\\1\end{pmatrix}: \quad r\in\mathbb{R}\right\} 
# $$
# 
# Since (geometric multiplicity) $1 < 2$ (algebraic multiplicity) we cannot diagonalise the matrix $A$. 
# 
# We can, however, find a matrix $P$ such that
# 
# $$
#  A = PJP^{-1}
# $$
# 
# with 
# 
# $$
#  J = \begin{pmatrix} ~2&~1\\0&2\end{pmatrix}.
# $$
# 
# For example, the matrix
# 
# $$
#  P = \begin{pmatrix} -1 & ~0\\1 & 1\end{pmatrix}
# $$
# 
# is such that $\displaystyle P^{-1}AP = \begin{pmatrix} ~2&~1\\0&2\end{pmatrix}$. 
# 
# To see this note that
# 
# $$
#  P^{-1} = \begin{pmatrix} -1&~0\\1&1\end{pmatrix}
# $$
# 
# and so 
# 
# $$
# P^{-1}AP = \begin{pmatrix} -1&~0\\1&1\end{pmatrix}\begin{pmatrix}1&-1\\1&3\end{pmatrix}\begin{pmatrix} -1 & ~0\\1 & 1\end{pmatrix} = \begin{pmatrix} ~2&~1\\0&2\end{pmatrix}.
# $$
# 
# More generally, it can be shown that given any real matrix $\displaystyle A\in\mathbb{R}^{n\times n}$, one can always decompose it as 
# 
# $$
#  A = PJP^{-1},
# $$
# 
# where
# 
# $$
# J = \begin{pmatrix} J_{n_1}(\lambda_1)&&\\&\ddots&\\&&J_{n_k}(\lambda_k)\end{pmatrix}, \qquad n_1+n_2+\cdots +n_k = n,
# $$ (JNF)
# 
# is a *block diagonal matrix* with blocks $\displaystyle J_{n_i}(\lambda_i)$ taking the form
# 
# $$
# J_{n_i}(\lambda_i) = \begin{pmatrix} \lambda_i&1&\\&\ddots&1\\&&\lambda_i\end{pmatrix} \in\mathbb{R}^{n_i\times n_i},
# $$
# 
# and $\displaystyle J_1(\lambda_i) = [\lambda_i]$. 
# 
# Notice that if each Jordan block $\displaystyle J_{n_i}(\lambda_i)$ is one-dimensional, that is, all $n_i=1$ and $k=n$, then the Jordan matrix $J$ is diagonal.
# 
# The matrix in {eq}`JNF` is the *Jordan normal (or canonical) form* of the matrix $A$.
# 
# <br>
# 
# **Example 2.0.2** Below we highlight how the above notation is used to distingiush between matrices that can be diagonalised and those that cannot.
# 
# $$
# J_1 = \begin{pmatrix} 2 &0 &0\\0 &2 &0\\0 &0 &2\end{pmatrix} = \begin{pmatrix} J_1(2)&&\\&J_1(2)&\\&&J_1(2) \end{pmatrix}, \quad J_2 = \begin{pmatrix} 2 &0 &0\\0 &2 &1\\0 &0 &2\end{pmatrix} = \begin{pmatrix} J_1(2)&\\&J_2(2) \end{pmatrix} 
# $$
# 
# and
# 
# $$
# J_3 = \begin{pmatrix} 2 &1 &0\\0 &2 &1\\0 &0 &2\end{pmatrix} = \left[J_3(2)\right].
# $$
# 
# The number of Jordan blocks corresponding to a given eigenvalue is the geometric multiplicity of the eigenvalue, which is the dimension of the associated eigenspace. 
# 
# It follows that the GMs for the matrices $J_1, J_2$ and $J_3$ are $3, 2$ and $1$, respectively.

# In[ ]:




