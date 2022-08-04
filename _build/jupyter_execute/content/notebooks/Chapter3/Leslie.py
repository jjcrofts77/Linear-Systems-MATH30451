#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/jjcrofts77/Linear-Systems-MATH30451/blob/main/content/notebooks/Chapter3/Leslie.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# In[ ]:





# # 3.4 Leslie Population Growth Models (optional)
# 
# In this section we consider the *Leslie model*, which was originally used to
# describe the growth of he female side of an animal population which could be
# divided into various classes, each of which was of equal duration. Then, if the
# maximum age of the population was $L$ and there were $n$ of these classes, the
# duration of each would be $L/n$ and it would be natural to observe the state of
# the populations at times $\displaystyle T_0=0, T_1=L/n, T_2=2L/n,\ldots,
# T_k=\frac{kL}{n},\ldots$
# 
# Suppose then that $\displaystyle b_i~(i=1,2,\ldots,n)$ is the average number of females born to animals in class $i$ and $\displaystyle s_i~(i=1,2,\ldots, n-1)$ is the fraction of class $i$ that survives to class $i+1$ (females in the $n$th class do not survive). Then if $x^{(i)}_k$ is the number of females in the class $i$ at the $k$th stage we see that, if 
# 
# $$
# \mathbf{x}_k=[x^{(1)}_k, x^{(2)}_k,\ldots,x^{(n)}_k]^T,
# $$
# 
# then
# 
# $$
# \mathbf{x}_{k+1} = L\mathbf{x}_k,
# $$
# 
# where
# 
# $$
# L=\begin{pmatrix} b_1&b_2&\ldots&b_{n-1}&b_n\\
# s_1&0&\ldots&0&0\\
# 0&s_2&\ldots&0&0\\
# \vdots&\vdots&&\vdots&\vdots\\
# 0&0&\ldots&s_{n-1}&0\end{pmatrix},
# $$
# 
# <br>
# 
# is a *Leslie matrix* which describes the evolution of the population.
# 
# For the above system we know that $\displaystyle\mathbf{x}_k=L^k\mathbf{x}_0$ and that the
# key to its long-term behaviour lies in the eigenvalues of the matrix $L$, which
# are given by $\displaystyle |L-\lambda I_n|=0$, where
# 
# $$
#  |L-\lambda I_n| &= \left|\begin{array}{cccccc}
#  b_1-\lambda&b_2&b_3&\ldots&b_{n-1}&b_n\\
#  s_1&-\lambda&0&\ldots&0&0\\
#  0&s_2&-\lambda&\ldots&0&0\\
#  \vdots&\vdots&\vdots&&\vdots&\vdots\\
#  0&0&0&\ldots&s_{n-1}&-\lambda
#  \end{array}\right|,\\
# &=(b_1-\lambda)\left|\begin{array}{ccccc}-\lambda&0&\ldots&0&0\\
# s_2&-\lambda&\ldots&0&0\\
# \vdots&\vdots&&\vdots&\vdots\\
# 0&0&\ldots&s_{n-1}&-\lambda\end{array}\right|-b_2
# \left|\begin{array}{ccccc}s_1&0&\ldots&0&0\\0&-\lambda&\ldots&0&0\\
# \vdots&\vdots&&\vdots&\vdots\\
# 0&0&\ldots&s_{n-1}&-\lambda\end{array}\right|+\cdots\\
# &+(-1)^{n-2}b_{n-1}\left|\begin{array}{ccccc}s_1&-\lambda&0&\ldots&0\\
# 0&s_2&-\lambda&\ldots&0\\\vdots&\vdots&\vdots&&\vdots\\
# 0&0&0&\ldots&-\lambda\end{array}\right|+
# (-1)^{n-1}b_n\left|\begin{array}{ccccc}s_1&-\lambda&0&\ldots&0\\
# 0&s_2&-\lambda&\ldots&0\\\vdots&\vdots&\vdots&&\vdots\\
# 0&0&0&\ldots&s_{n-1}\end{array}\right|,\\\nonumber
# &=(b_1-\lambda)(-\lambda)^{n-1}-b_2s_1(-\lambda)^{n-2}
# +\cdots+(-1)^{n-2}b_{n-1}s_1s_2\cdots s_{n-2}(-\lambda)+\cdots\\
# &+(-\lambda)^{n-1}b_ns_1s_2\cdots s_{n-1},\\
# &=
# (-1)^n\left[\lambda^n-b_1\lambda^{n-1}-b_2s_1\lambda^{n-2}-\cdots - b_{n-1}s_1s_2\cdots s_{n-2}\lambda-b_ns_1s_2\cdots s_{n-1}\right],\\
# &=(-1)^n\lambda^n\left[1-q(\lambda)\right],
# $$
# 
# where
# 
# $$
# q(\lambda) =
# \frac{b_1}{\lambda}+\frac{b_2s_1}{\lambda^2}+\frac{b_3s_1s_2}{\lambda^3}
# +\cdots+\frac{b_{n-1}s_1s_2\cdots s_{n-1}}{\lambda^{n-1}}+\frac{b_ns_1s_2\cdots
# s_{n-1}}{\lambda^n}.
# $$
# 
# Then, eigenvalues of $L$ are given by $q(\lambda)=1$ and, as the graph of
# $q(\lambda)$ shows, there is only one, unrepeated, positive eigenvalue $\lambda$. Eigenvectors corresponding to this eigenvalue have the form
# $\displaystyle\mathbf{u}=[u_1,u_2,\ldots,u_n]^T$, where
# 
# $$
# (L-\lambda_1I_n)\mathbf{u}=\begin{pmatrix} b_1-\lambda_1&b_2&b_3&\ldots&b_{n-1}&b_n\\
# s_1&-\lambda_1&0&\ldots&0&0\\
# 0&s_2&-\lambda_1&\ldots&0&0\\
# \vdots&\vdots&&\vdots&\vdots\\
# 0&0&0&\ldots&s_{n-1}&-\lambda_1\end{pmatrix}\begin{pmatrix}
# u_1\\u_2\\u_3\\\vdots\\u_{n-1}\\u_n\end{pmatrix}=\mathbf{0}.
# $$
# 
# <br>
# 
# **Example 3.4.1** Let 
# 
# $$
# L_1=\begin{pmatrix}
# 0&7&8\\\frac{1}{2}&0&0\\0&\frac{1}{4}&0\end{pmatrix},
# $$
# 
# then
# 
# $$
# |L_1-\lambda I_3| = \left|\begin{matrix}
# -\lambda&7&8\\\frac{1}{2}&-\lambda&0\\0&\frac{1}{4}&-\lambda
# \end{matrix}\right|=-\lambda^3+\frac{7\lambda}{2}
# +1=-(\lambda-2)(\lambda^2+2\lambda+\frac{1}{2}).
# $$
# 
# Eigenvalues of $L$ are, therefore, $2$ and $-1\pm\frac{1}{\sqrt{2}}$.
# 
# Eigenvectors corresponding to the eigenvalue $2$ are given by
# 
# $$
# (L-2I_3)\mathbf{u} = 0, \quad\text{i.e.,}~
# \begin{pmatrix}-2&7&8\\\frac{1}{2}&-2&0\\0&\frac{1}{4}&-2\end{pmatrix}
# \begin{pmatrix} u_1\\u_2\\u_3\end{pmatrix}=0.
# $$
# 
# Solving gives the eigenvector $[4,1,0.125]$. We have thus shown that a population which occurs in the ratio $32:8:1$ will double steadily.
# 
# Now, in the example above the single positive eigenvalue is, in fact, strictly
# dominant. Generally, this eigenvalue is always a dominant one but it is not
# necessarily the only such eigenvalue.
# 
# <br>
# 
# **Example 3.4.2** For 
# 
# $$
# L_2=\begin{pmatrix}
# 0&0&6\\\frac{1}{2}&0&0\\0&\frac{1}{3}&0\end{pmatrix},
# $$
# 
# we have
# 
# $$
# |L_2-\lambda I_3| = \left|\begin{matrix}
# -\lambda&0&6\\\frac{1}{2}&-\lambda&0\\0&\frac{1}{3}&-\lambda
# \end{matrix}\right|=-\lambda^3+1=-(\lambda-1)(\lambda^2+2\lambda+1).
# $$
# 
# This matrix has eigenvalues $\displaystyle 1, (-1\pm i\sqrt{3})/2$ and these all have magnitude 1.
# 
# Note that in this case the Cayley-Hamilton Theorem tells us that $\displaystyle L^3_2=I$, so that we can see that $\displaystyle\mathbf{x}_{k+3}=L^3\mathbf{x}_k=\mathbf{x}_k$ for any $k$, *i.e.*, in this situation the population evolves periodically, with the same situation occurring at intervals of $3$ times the duration of the classes.
# 
# Interestingly enough, there is a condition which can be used to test if the
# single positive eigenvalue of a given Leslie matrix $L$ is strictly dominant and
# this is simply that it should have two successive elements in its first row
# which are non-zero. In other words, there should be two successive fertile
# classes in the situation previously described. Note that, as we should expect,
# the matrix $L_1$ in Example 3.4.1 satisfies this condition but $L_2$ in Example
# 3.4.2 does not.
# 
# When a Leslie matrix satisfies the condition which guarantees that it has a 
# strictly dominant positive eigenvalue, this eigenvalue can be found by the power
# method, which ensures that from a general starting point we eventually reach a 
# situation in which $\displaystyle\mathbf{x}_{k+1}=\lambda_i\mathbf{x}_k$. (i.e., eventually the population reaches a state where the number in each class changes steadily at a rate of $\displaystyle 100(\lambda_i-1)\%$.)
