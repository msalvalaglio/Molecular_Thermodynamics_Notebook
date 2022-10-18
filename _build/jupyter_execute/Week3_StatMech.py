#!/usr/bin/env python
# coding: utf-8

# # An Introduction to Statistical Mechanics
# 
# __Defining statistical mechanics:__
# Statistical Mechanics provides a connection between microscopic motion of individual atoms of matter and macroscopically observable properties such as temperature, pressure, entropy, free energy, heat capacity, chemical potential, viscosity, spectra, reaction rates, etc.
# 
# __Why do we need Statistical Mechanics?__
# 
# - Statistical Mechanics provides the microscopic/fundamental basis for thermodynamics, which, otherwise is a phenomenological theory.
# 
# - Microscopic details allows access to properties that are typically not accessible from classical thermodynamics - a typical example is the microscopic stucture of crystal nuclei, proteins, micelles etc... .
# 
# - Statistical mechanics is based on the detailed microscopic description of a system. Statiustical mechanics allows to connect models of the microscopic world of atoms/particles to macroscopic thermodynamic quantities, enabling the interpretation and elucidation of experimental results. 
# 
# 
# __Do we need atomic detail?__ (Yes and No).
# 
# _Yes_ – There are plenty of physical-chemical phenomena and processes for which we really would like to know how individual atoms are moving as a process evolves. From a theoretical point of view, although we cannot follow the behavior of every particle in a macroscopic system ($\approx{10^{23}}$ particles), we can model the motion of every particle in a system containing $10^3-10^6$ particles - which might capture most of the features of true macroscopic matter. 
# 
# _No_ – Intuitively, we would expect that if we were to follow the evolution of a large number of systems all described by the same set of forces but having starting from different initial conditions, these systems would have essentially the same macroscopic characteristics, e.g. the same temperature, pressure, etc. even if the microscopic detailed evolution of each system in time would be very different. This idea suggests that the microscopic details are largely unimportant.
# However, microscopic details are largely unimportant when discussing and obtaining macroscopic properties. This suggests that macroscopic quanities can be obtained by averging out microscopic characteristics via averaging. This is idea was developed by the founders of Thermodynamics and Statistical Mechanics:  Gibbs, Maxwell, and Boltzmann.
# 
# __Ensemble:__ An ensemble consists of large number of systems, each described by the same set of microscopic coordinates and forces and sharing some common macroscopic property (e.g. their total energy). Each system is assumed to evolve under the microscopic laws of motion from a different set of initial conditions so each system will evolve differently with respect from all the others.  
# 
# __Ensemble Average:__ At the heart of the ensemble concept there is the idea of ensemble average. It states that macroscopic observables can be calculated by averaging properties over the systems in an ensemble. 
# The ensemble averaging process can be translated mathematically as follows. Let's denote with $A$ a macroscopic, time-independent, property and with $a$ a microscopic function that used to compute $A$. As an example, $A$ could be the temperature, and $a$ the kinetic energy (a microscopic function of the particle velocities). The ensemble average definition allows to define $A$ as the average of $a$ computed from all systems belonging to the ensemble:
# 
# $$
# A=\frac{1}{N}\sum_{i=1}^Na_i
# $$
# 
# where $N$ is the total number of members of the ensemble and $a_i$ is the value of a in the $i^{th}$ system: 
# 
# The questions that naturally arise are: 
# 
# - How do we construct an ensemble?
# 
# - How do we perform averages over an ensemble?
# 

# __Phase space, Hamiltonian__ 
# 
# A microscopic configuration of a system of N particle is fully determined by its 3N coordinates $\mathbf{x}$ and 3N momenta $\mathbf{p}$. 
# Phase space is, the 6N dimensional Cartesian space where each of the 6N coordinates and momenta is assigned to one of 6N mutually orthogonal axes. 
# A point in phase space is specified by a specific realisation of the 6N coordinates and momenta, and can be indicated by the 6N dimensional vector: 
# 
# $$
# \mathbf{x} = (p1,...,pN,r1,...,rN)
# $$
# 
# The time evolution or trajectory of a system can be expressed as $x(t)$: the time evolution of the phase space vector.
# 
# The total energy of a microscopic system corresponds to its Hamiltonian, (i.e. the sum of the Kinetic and Potential energies).
# The Hamiltonian of a system is a function of its corresponding phase space vector: $H(\mathbf{x}(t))$
# 
# The law of Energy conservation can thus be expressed as a condition on the phase space vector: 
# 
# $$
# H(\mathbf{x}(t)) = const = E
# $$
# 
# This constraint implies that the time evolution of a system at constant total energy takes place in a 6N − 1 dimensional hypersurface subject to this constraint.  
# 
# 
# __Classical Microstate__
# 
# A microscopic state, or microstate, of a classical system is a complete specification of the 6N coordinates that define $\mathbf{x}$. At constant Energy any valid microstate lies on the constant energy hypersurface, $H(\mathbf{x}) = E$. 
# The concept of classical microstates allows us to elaborate a more formal definition of an __ensemble__.
# An _ensemble is a collection of systems sharing one or more macroscopic characteristics but each being in a unique microstate_. The complete ensemble is specified by all microstates consistent with the common macroscopic characteristics of the ensemble (i.e. constant Energy).
# 
