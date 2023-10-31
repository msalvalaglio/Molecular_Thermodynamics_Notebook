#!/usr/bin/env python
# coding: utf-8

# # The Microcanonical Ensemble
# 
# The microcanonical ensemble is built upon the so called postulate of equal a priori probabilities:
# 
# _Postulate of equal a priori probabilities: For an isolated macroscopic system in equilibrium, all microscopic states corresponding to the same set of macroscopic observables are equally probable._
# 
# To intuitively absorb the idea of ensemble let's consider a system and its phase space vector $\mathbf{x}(t)$, evolving in time at constant total energy, volume and number of particles. The configurations sampled, i.e. the realizations of $\mathbf{x}$ that constitute $\mathbf{x}(t)$ are sampling the the constant energy hypersurface and are therefore equally probable. 
# 
# The number of equally probable states $\Omega(N,V,E)$ defines the __partition function__ for the Microcanonical ensemble. 
# 
# The partition function of the microcanonical ensemble is proportional to the integral in phase space over all the states at constant Energy: 
# 
# $$
# \Omega(N,V,E)\propto\int{d\mathbf{x}}\delta(H(\mathbf{x})-E)
# $$
# 
# the proportionality constant is equal to: $C_N=E_0/N!h^{3N}$
# 
# Here $h$ is a constant with units Energy$\times$Time, and  $E_0$ is a constant having units of energy. The extra factor of  $E_0$ is needed because the $\delta$ function has units of inverse energy. Such a constant however has no effect at all on any properties). Thus,  $\Omega(N,V,E)$ is non-dimensional. 
# 
# The microcanonical partition function quantify the number of microstates accessible to a system at constant energy. 
# Boltzmann identified $\Omega(N,V,E)$ as the microscopic quantity defining $S$ - the Entropy - which we was introduced as a postulate while stating the second law of Thermodynamics. As $S$,  $\Omega(N,V,E)$ is a natural function of $N$, $V$ and $E$. 
# 
# The famous Boltzmann's relation between $S$ and  $\Omega(N,V,E)$ is: 
# 
# $$
# S(N,V,E)=k\ln{\Omega(N,V,E)}
# $$
# 
# where k is Boltzmann's constant. 
# 
# This is a key relation, as it establishes a connection between the macroscopically measurable thermodynamic properties of a system and its microscopic details.
# 
# In particular, by recalling the relationships introduced in Week 1, defining T and the forces associated with the extensive mechanical variables we get: 
# 
# $$
# k\left(\frac{\partial{\ln{\Omega}}}{\partial{V}}\right)_U=\left(\frac{\partial{S}}{\partial{V}}\right)_U=-\left(\frac{\partial{S}}{\partial{U}}\right)_V=\frac{P}{T}
# $$
# 
# $$
# k\left(\frac{\partial{\ln{\Omega}}}{\partial{n_i}}\right)_U=\left(\frac{\partial{S}}{\partial{n_i}}\right)_U=-\left(\frac{\partial{S}}{\partial{U}}\right)_{n_i}=-\frac{\mu_i}{T}
# $$
# 
# and
# 
# $$
# k\left(\frac{\partial{\ln{\Omega}}}{\partial{U}}\right)_{n_i,V}=\left(\frac{\partial{S}}{\partial{U}}\right)_{n_i,V}=\frac{1}{T}
# $$
# 
