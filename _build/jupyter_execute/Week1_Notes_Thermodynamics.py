#!/usr/bin/env python
# coding: utf-8

# # Thermodynamics Revision, with an eye to Statistical Mechanics 
# 
# ## Recommended reading  
# - Chapters 1, 2 Prausnitz 
# - Chapters 1, 2 Chandler (Introduction to Modern Statistical Mechanics)
# 
# ## Objective of these notes: 
# To guide revision of the principles of Thermodynamics to:
# 
# - Revise the first and second law of thermodynamics
# - Define the concept of thermodynamic equilibrium 
# - Define the concept of chemical potential
# - Define state functions and thermodynamic potentials
# 
# 

# ## First Law of Thermodynamics. 
# 
# The first law of thermodynamics is about internal __Energy__. It is based on two postulates: 
# 1. Internal Energy is an __extensive property__.
# 2. Internal Energy is __conserved__.
# 
# This leads to the formulation of the first law of thermodynamics: 
# 
# $$
# dU=\mathbf{d}Q+\mathbf{d}W
# $$
# 
# where $Q$ and $W$ represent heat and work, and $\mathbf{d}$ indicate inexact differentials. 
# 
# $\mathbf{d}W$ is the differential work done by manipulating mechanical constraints. It can be written generally as: 
# 
# $$
# dW=\mathbf{f}{d}\mathbf{X}
# $$
# 
# where $\mathbf{f}$ represents a vector of applied forces and $\mathbf{X}$ a vector of their corresponding mechanical extensive variables. 
# 
# It should be noted that heat and work only represent forms of energy transfer, and once the tranfer has taken place there is no $Q$ or $W$ quantity remains associated with a given state. This is at why $\mathbf{d}Q$ and $\mathbf{d}W$ are _inexact_ differentials.
# 

# ## Second Law of Thermodynamics
# 
# The second law of thermodynamics concerns with the reversibility of transformations, and emerges from assmptions on the nature of equlibrium states. 
# Following Chandler's approach we start from a postulate of the second law: 
# 
# - There is an extensive function of state, $S(U, X)$, which is a monotonically increasing function of U, and if state B is adiabatically accessible from state A, then $S_B\geq{S}_A$.
# 
# This implies that the change $\Delta{S} = S_B - S_A$ is zero for a reversible adiabatic ($dQ=0$) process, and otherwise $\Delta{S}$ is positive for any natural irreversible adiabatic process.
# 
# The extensive function of state, $S(U, X)$, is the entropy of the system. Hence, the entropy change for a reversible adiabatic process is null. Note also that entropy is a function of state. That means it is defined for those states characterized by $U$ and $\mathbf{X}$.
# 
# ### Definition of Temperature: 
# Since $S$ is a monotonically increasing function of $U$, $\rightarrow$ $\left(\frac{\partial{S}}{\partial{U}}\right)_X>{0}$, or $\left(\frac{\partial{U}}{\partial{S}}\right)_X\geq{0}$. This derivative defines the absolute temperature $T$: 
# 
# $$
# \left(\frac{\partial{U}}{\partial{S}}\right)_X\equiv{T}
# $$

# ### Mathematical formulation of the Second Law: the first step towards deriving equilibrium criteria for processes at constant Temperature  
# Let's start by writing the differential of the Entropy function, introduced above: 
# 
# $$
# dS(U,X)=\left(\frac{\partial{S}}{\partial{U}}\right)_XdU + \left(\frac{\partial{S}}{\partial{X}}\right)_UdX
# $$
# 
# This indicates a variaton of entropy associated to any variation of internal energy $dU$, and of the vector of extensive mechanical properties $X$.
# For any reversible transition $dU=\mathbf{d}Q_{rev}+f\mathbf{d}X$, and $dS=0$, hence: 
# 
# $$
# 0=\left(\frac{\partial{S}}{\partial{U}}\right)_X(\mathbf{d}Q_{rev}+f\mathbf{d}X)+\left(\frac{\partial{S}}{\partial{X}}\right)_UdX
# $$
# 
# This equation is true for any reversible transition, including the adiabatic ones. Hence it must be that: 
# 
# $$
# 0=\left[\left(\frac{\partial{S}}{\partial{U}}\right)_Xf+\left(\frac{\partial{S}}{\partial{X}}\right)_U\right]\mathbf{d}X
# $$
# 
# and 
# 
# $$
# \left(\frac{\partial{S}}{\partial{X}}\right)_U=-\left(\frac{\partial{S}}{\partial{U}}\right)_Xf
# $$
# 
# By introducing the definition of absolute temperature gives: 
# 
# $$
# \frac{f}{T}=\left(\frac{\partial{S}}{\partial{X}}\right)_U
# $$
# 
# This expression enables to rewrite the differential of entropy $dS$ as: 
# 
# $$
# dS=\frac{1}{T}dU - \frac{f}{T}dX
# $$
# 
# which, for a pure homogeneous phase leads to the fundamental expression:
# 
# $$
# dU=TdS-PdV
# $$
# 
# _N.B._ This equation establishes the triplet of state variables $U$, $S$ and $V$ as a so-called fundamental grouping. Different choices are possible, leading to the definition of different fundamental groupings with their own thermodynamic potential. 
# 

# ## Changing the fundamental grouping from V, S, U: Legendre transform and Thermodynamic potentials
# 
# So far we have discussed of Energy, Entropy and Volume - and learned that the equilibrium conditions for a system as a function of S and V correspond to minimum E. 
# If one needs to characterise equilibrium as a function of a different set of variables than entropy and volume, can a therodynamic function that plays the same role of Energy be identified? 
# 
# The answer is, yes. Such function, hereafter indicated as thermodynamic potential, can be obtained by swapping variables by performing a _Legendre Transform_. 
# 
# The simplest way of demonstrating the process of operating a Legendre transform to obtain a function of a new set of variables is to start from original function, and subtract the product of the variables that need to be swapped. Differentiating such new function allows to derive a new fundamental relationship defining a thermodynamic potential for the new grouping. 
# 
# ### Swapping T and S, from Energy to Helmholtz Free Energy. 
# 
# If we define a new grouping in which the entropy and temperature are swapped, i.e. we want to define a new thermodynamic potential able to characterise equilibrium for a system at constant T and V we should start by defining a new function as follows: 
# 
# $$
# A=U-TS
# $$
# 
# where $A$ is the new thermodynamic potential, the Helmholtz free energy, $T$ is the temperature and $S$ the entropy. Computing the differential for A give us: 
# 
# $$
# dA=dU-SdT-TdS=TdS-PdV-SdT-TdS=-SdT-PdV
# $$
# 
# which is the new thermodynamic relationship that allows to characterise equilibrium as a function of $V$ and $T$, rather than $V$ and $S$.  
# 
# ### Swapping P and V, from Energy to Enthalpy. 
# 
# Let's start by Legendre transforming the $U$ to swap $V$ with its intensive conjugate variable $P$, this yields H, or enthalpy: 
# 
# $$
# H=U+PV
# $$
# 
# the fundamental equation, in a system at constant P and S then becomes: 
# 
# $$
# dH=dU+VdP+PdV=TdS-PdV+VdP+PdV=TdS+VdP
# $$
# 
# ### Gibbs Free Energy
# 
# Gibbs Free Energy is obtained by swapping both T/S and P/V. It is defined as: 
# 
# $$
# G=U-TS+PV
# $$
# 
# the expression of the Gibbs free energy differential $dG$ is thus: 
# 
# $$
# dG=dH-SdT-TdS=TdS+VdP-SdT-TdS=-SdT+VdP
# $$
# 

# # The Phase Equilibrium Problem
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
# 
# $$
# S=S^\alpha+S^\beta
# $$
# 
# volume 
# 
# $$
# V=V^\alpha+V^\beta
# $$
# 
# and the number of moles of each i$^{th}$ component in the different phases: 
# 
# $$
# n_i=n_i^\alpha+n_i^\beta
# $$
# 
# Consider now to define the equilibrium conditions for a system where the total $S$, $V$ and $n_i$ are constant, but these quantities can be redistributed between phases $\alpha$ and $\beta$. This requires that: 
# 
# $$
# \delta{S}^{\alpha} = - \delta{S}^{\beta} 
# $$
# 
# $$
# \delta{V}^{\alpha} = - \delta{V}^{\beta} 
# $$
# 
# $$
# \delta{n_i}^{\alpha} = - \delta{n_i}^{\beta} 
# $$
# 
# The equilibrium can be investigated with a variational principle, i.e. by introducing a perturbation $\delta{U}$ from equilibrium. If the system is in equilibrium the perturbation can only be associated to and _increase_ in U: 
# 
# $$
# \delta{U}_{S,V,n_i}\geq 0 = (T^\alpha-T^\beta)\delta{S^{\alpha}}-(p^{\alpha}-p^{\beta})\delta{V^{\alpha}}+(\mu_i^{\alpha}-\mu_i^{\beta})\delta{n_i^{\alpha}}
# $$
# 
# Since the variations in $n$, $S$, and $V$ are all decoupled, and can be both positive or negative, the only valid solution for the equilibrium condition stated above is that the temperature, pressure, and chemical potential of component i remain constant across the phase boundary between $\alpha$ and $\beta$, namely:
# 
# $$
# {T}^{\alpha} = {T}^{\beta} 
# $$
# 
# $$
# {p}^{\alpha} = {p}^{\beta} 
# $$
# 
# $$
# {\mu_i}^{\alpha} = {\mu_i}^{\beta} 
# $$
# 
# ## Degrees of freedom of a single multicomponent phase.
# 
# Given what we have seen so far, can we determine how many degrees of freedom are associated to a single phase?  
# 
# In a single homogeneous multicomponent phase at constant $T$ and $P$, containing $m$ component the number of degrees of freedom is: 
# 
# $$
# DOF=m+2-c
# $$
# 
# where $c$ is the number of constraints. For a single phase there is only one constrain that applies, the Gibbs Duhem equation. Hence $DOF=m+1$.
# 
# ### The Gibbs Duhem Equation
# The fundamental equation emerging from the first and second laws of Thermodynamics for a multicomponent system reads:
# 
# $$
# dU=TdS-pdV+\sum_i \mu_idn_i
# $$
# 
# Integrating this equation to obtain $U$ along from a state with no molecules to a state where a final number $n_i$ of molecules is present, along an isothermal, isobaric transformation at constant compositionleads to the following expression for $U$:
# 
# $$
# U=TS-pV+\sum_i\mu_in_i
# $$
# 
# Differentiating $U$ leads to: 
# 
# $$
# dU=TdS+SdT-pdV-VdP+\sum_i\mu_idn_i+\sum_id\mu_in_i
# $$
# 
# Equating the right hand sides of these two expressions leads to: 
# 
# $$
# TdS-pdV+\sum_i \mu_idn_i=TdS+SdT-pdV-VdP+\sum_i\mu_idn_i+\sum_id\mu_in_i
# $$
# 
# which implies that the following equality is always satisfied for an homogenous multicomponent phase: 
# 
# $$
# SdT-Vdp+\sum_{i}n_{i}d\mu_i=0
# $$
# 
# This equality is the Gibbs-Duhem equation, which allows express the chemical potential as a function of the temperature and pressure. 
# 
# 
# ## Degrees of fredom in a multiphase, multicomponent system: The Gibbs phase rule
# 
# Let's consider a system comprising $\nu$ homogenous phases at equilibrium, each containing $m$ components. The number of degrees of freedom associated to each phase is: 
# 
# $$
# DOF=\nu(m+1)
# $$
# 
# However at equlibrium we have that for every couple of phases ($\nu-1$), $T$, $P$ and all $m$ chemical potentials needs to be constant across phase boundaries. 
# This means that we have a number of constraints equal to: 
# 
# $$
# c=(\nu-1)(m+2)
# $$
# 
# This leads to a number of degrees of freedom equal to: 
# 
# $$
# DOF=\nu(m+1)-(\nu-1)(m+2)=m+2-\nu
# $$
