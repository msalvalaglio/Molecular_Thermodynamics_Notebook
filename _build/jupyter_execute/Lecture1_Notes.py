#!/usr/bin/env python
# coding: utf-8

# # Topic 2: The Phase Equilibrium Problem
# 
# ## References: 
# - Chapters 1, 2 Prausnitz 
# - Chapters 1, 2 Chandler
# 
# ## The problem: 
# _Relate quantitatively the variables that describe the state of equilibrium between two or more homogeneous phases that are free to exchange energy and matter._
# 
# ## Definitions: 
# - __Homogeneous phase at equilibrium__: Region of space where the intensive properties are constant. 
# - __Intensive properties__: properties that are independent of mass, size, or shape of the phase (i.e. T, P, r, x, y, etc.)

# ## Equilirbium in a multiphase, multicomponent, system: 
# 
# Consider a system where: 
# - Two or more phases are present
# - Each phase can be considered an open system within the over all closed system
# 
# How do we know when such a system is at equilibrium?
# 
# The energy, $U$, is an extensive property of the system. In an heterogeneous multipase system containing at least two phases $\alpha$ and $\beta$ the total energy of the system is $U=U^\alpha+U^\beta$. The same can be written for the other extensive variables defining the state of the multiphase system such as entropy: 
# \begin{equation}
# S=S^\alpha+S^\beta
# \end{equation}
# volume 
# \begin{equation}
# V=V^\alpha+V^\beta
# \end{equation}
# and the number of moles of each i$^{th}$ component in the different phases: 
# \begin{equation}
# n_i=n_i^\alpha+n_i^\beta
# \end{equation}
# 
# 
# Consider now to define the equilibrium conditions for a system where the total $S$, $V$ and $n_i$ are constant, but these quantities can be redistributed between phases $\alpha$ and $\beta$. This requires that: 
# \begin{equation}
# \delta{S}^{\alpha} = - \delta{S}^{\beta} 
# \end{equation}
# \begin{equation}
# \delta{V}^{\alpha} = - \delta{V}^{\beta} 
# \end{equation}
# \begin{equation}
# \delta{n_i}^{\alpha} = - \delta{n_i}^{\beta} 
# \end{equation}
# 
# The equilibrium can be investigated with a variational principle, i.e. by introducing a perturbation $\delta{U}$ from equilibrium. If the system is in equilibrium the perturbation can only be associated to and _increase_ in U: 
# 
# \begin{equation}
# \delta{U}_{S,V,n_i}\geq 0 = (T^\alpha-T^\beta)\delta{S^{\alpha}}-(p^{\alpha}-p^{\beta})\delta{V^{\alpha}}+(\mu_i^{\alpha}-\mu_i^{\beta})\delta{n_i^{\alpha}}
# \end{equation}
# 
# Since the variations in $n$, $S$, and $V$ are all decoupled, and can be both positive or negative, the only valid solution for the equilibrium condition stated above is that the temperature, pressure, and chemical potential of component i remain constant across the phase boundary between $\alpha$ and $\beta$, namely:
# \begin{equation}
# {T}^{\alpha} = {T}^{\beta} 
# \end{equation}
# \begin{equation}
# {p}^{\alpha} = {p}^{\beta} 
# \end{equation}
# \begin{equation}
# {\mu_i}^{\alpha} = {\mu_i}^{\beta} 
# \end{equation}
# 
# 
# ## Degrees of freedom of a single multicomponent phase.
# 
# Given what we have seen so far, can we determine how many degrees of freedom are associated to a single phase?  
# 
# In a single homogeneous multicomponent phase at constant $T$ and $P$, containing $m$ component the number of degrees of freedom is: 
# \begin{equation}
# DOF=m+2-c
# \end{equation}
# 
# where $c$ is the number of constraints. For a single phase there is only one constrain that applies, the Gibbs Duhem equation. Hence $DOF=m+1$.
# 
# ### In Class Activity 5: 
# Derive the Gibbs-Duhem equation: 
# \begin{equation}
# SdT-Vdp+\sum_{i=1}^{n_{species}}n_i\mu_i=0
# \end{equation}
# 
# ## Degrees of fredom in a multiphase, multicomponent system: The Gibbs phase rule
# 
# Let's consider a system comprising $\nu$ homogenous phases at equilibrium, each containing $m$ components. The number of degrees of freedom associated to each phase is: 
# \begin{equation}
# DOF=\nu(m+1)
# \end{equation}
# 
# However at equlibrium we have that for every couple of phases ($\nu-1$), $T$, $P$ and all $m$ chemical potentials needs to be constant across phase boundaries. 
# This means that we have a number of constraints equal to: 
# \begin{equation}
# c=(\nu-1)(m+2)
# \end{equation}
# 
# This leads to a number of degrees of freedom equal to: 
# \begin{equation}
# DOF=\nu(m+1)-(\nu-1)(m+2)=m+2-\nu
# \end{equation}

# ### Study Questions
# 
# - What do gradients in chemical potential cause when two phases are at contact with each other?
# - A system containing ethanol and water is split into liquid and vapor phases; the chemical potential of water in the liquid phase is different than the chemical potential of ethanol in the vapor phase: is the system at equilibrium?
