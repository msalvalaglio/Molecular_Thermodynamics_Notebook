#!/usr/bin/env python
# coding: utf-8

# # The Canonical Ensemble
# 
# In the canonical ensemble, all microstates are characterised by the same $N$, $V$ and $T$. To obtain the relevant thermodynamic potential for the canonical ensemble we can start from the microcanonical one. 
# 
# In the microcanonical ensemble, the entropy $S$ is a natural function of $N$,$V$ and $E$, i.e., $S=S(N,V,E)$. 
# This can be inverted to give the energy $E$ as a function of $N$,$V$, and $S$, i.e., $E=E(N,V,S)$. 
# As discussed in Week 1, by using the Legendre transformation to swap $S$ with $T$ one can obtain a new potential, called the Hemlholtz free energy, and is given the symbol $A(N,V,T)$. 
# 
# The Helmoholtz Free Energy is the fundamental energy in the canonical ensemble.

# __The Canonical Partition function__ 
# The question now becomes, what is the partition function in the canonical ensemble? 
# 
# The canonical partition function can be written as: 
# 
# $$
# Q(N,V,T) \propto \int{ d { \mathbf{x}}e^{-\frac{H(\mathbf{x})}{kT}}}
# $$
# 
# where $H(\mathbf{x})$ is the Hamiltonian of the system. 
# 
# The proportionality constant in the case of the Canonical partition function is $(N!h^{3N})^{-1}$ leading to the equality: 
# 
# $$
# Q(N,V,T) = \frac{1}{N!h^{3N}} \int{ d { \mathbf{x}}e^{-\frac{H(\mathbf{x})}{kT}}}
# $$
# 
# The fundamental relation between $A$ and $Q(N,V,T)$ is: 
# 
# $$
# A=kT\ln{Q(N,V,T)}
# $$
# 
# This equation provides a link between the microscopic and macroscopic variables at constant $N$, $V$ and $T$. 

# In[ ]:




