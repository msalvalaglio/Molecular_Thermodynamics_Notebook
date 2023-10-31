#!/usr/bin/env python
# coding: utf-8

# # From Intermolecular Interactions to Dynamics at the Atomistic Scale
# 
# ## Molecular Dynamics, Newton's Second Law and The Equations of Motion
# In MD simulations, atoms are treated as point particles evolving in time and space by following the principles of classical mechanics.  
# Newton's Second Law describes how the velocity of a particle changes when it is subjected to an external force, and in particular - when a force of magnitude $\textbf{F}$ is exerted on a body $i$, it results in acceleration $\textbf{a}$ according to:
# 
# $$
# 	\textbf{a}=\frac{\textbf{F}}{m_i}
# $$
# 
# where $m_i$ is the mass of the given particle.
# 
# At any given time $\textit{t}$ the position of A in Cartesian space is defined by the vector $\textbf{r}(t)=\left[x(t),y(t),z(t)\right]$. 
# 
# The velocity of particle $i$, $\textbf{v}(t)$, is the derivative of the position vector with respect to time, while the acceleration $\textbf{a}(t)$ is the derivative of the velocity with respect to time. 
# 
# For a system containint N-particles, the equations of motion read:
# $$
# \textbf{v}_i=\frac{\text{d}\textbf{r}_i}{\text{d}t}
# $$
# 
# $$
# \textbf{F}_i={m_i}{\frac{\text{d}\textbf{v}_i}{\text{d}t}}
# $$
# 
# Given a set of initial positions $\{\mathbf{r_1}(0),\mathbf{r}_2(0)...\mathbf{r}_N(0)\}$ and velocities $\{\mathbf{v_1}(0),\mathbf{v}_2(0)...\mathbf{v}_N(0)\}$, the time evolution of the system can be computed by solving the equations of motion. 
# 
# An exact solution of the equations cannot be computed in the general case of system composed on $N$ particles system, therefore numerical schemes, commonly based on Taylor series expansion are used. 

# ## MD Simulations Algorithm 
# 
# The MD algorithm consists of a few conceptual steps: 
# 
# - system initialisation
# - computing of the forces
# - integrating the equations of motion 
# - computation and output of ensemble properties.
# 
# ### System Initialization
# Initialising the system requres defining/identifying: 
# - A Force Field (FF). A FF is a mathematical formulation of the interaction potential $V$, complete with parameters that capture the physics of specific interatomic interactions. It includes functions accounting for non-bonded interactions, such as the Lennard-Jones potential, and bonded interaction functions, such as the harmonic potential. The force field is an important component of an MD simulation, and force fields are validated by comparing system properties to experimental data or \textit{ab-initio} calculations to ensure a physically meaningful system representation. 
# 
# 
# - The initial conditions: initial positions and velocities to all particles in the system.
# - Parameters such as number of particles, temperature, time step etc., as well as assigning   
# 
# 
# ### Forces Calculation  
# The second step of the algorithm involves calculating the force acting on every particle according to:
# 
# $$
# \textbf{F}_i = -\frac{\partial V}{\partial \textbf{r}_i} = \bf{F}_{i,bonded} + \sum_j{ \bf{F}_{ij,non-bonded}} + \bf{F}_{ext}
# $$
# 
# 
# The force is computed each time step $\Delta t$. Typically, for computational efficiency, the contribution to the force on particle $i$ due to all its neighbours is considered within a cut-off radius. 
# 
# The term ${F}_{i,bonded}$ refers to intramolecular interactions, i.e. interactions between atoms belonging to the same molecule and determine the intramolecular motions such as bond stretching, angle bending, and dihedral torsional rotations. As the atoms involved in a bond, angle or dihedral are considered to oscillate around the equilibrium position, these motions are usually described by a harmonic potential.
# 
# ### Dynamics Propagation
# After the forces have been obtained, the movement of the atoms is propagated by numerically solving the equations of motion. A number of algorithms have been designed to integrate the equations of motion which provide different level of accuracy depending on their application. 
# Commonly used are methods based on the Verlet algorithm. This algorithm is a combination of Taylor expansions, allowing to compute the atomic positions at time $(t+\Delta t)$ based on the atomic positions at time $t$ and $(t-\Delta t)$. The key equation of the Verlet algorithm is:  
# 
# $$
# r\left(t+\Delta t\right) = 2r(t) -r(t-\Delta{t}) +\frac{{\Delta{t}}^2}{m}\textbf {F}(t)
# $$
# 
# 
