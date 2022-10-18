#!/usr/bin/env python
# coding: utf-8

# In[1]:


__Canonical Ensemble__

Basic Thermodynamics

In the microcanonical ensemble, the entropy S is a natural function of N,V and E, i.e., S=S(N,V,E). This can be inverted to give the energy as a function of N,V, and S, i.e., E=E(N,V,S). Consider using Legendre transformation to change from S to T using the fact that

displaymath135

The Legendre transform  tex2html_wrap_inline606 of E(N,V,S) is

eqnarray140

The quantity  tex2html_wrap_inline610 is called the Hemlholtz free energy and is given the symbol A(N,V,T). It is the fundamental energy in the canonical ensemble.

The differential of A is

displaymath146

However, from A=E-TS, we have

displaymath154

From the first law, dE is given by

displaymath156

Thus,

displaymath158

Comparing the two expressions, we see that the thermodynamic relations are

eqnarray160


# The partition function
# 
# Consider two systems (1 and 2) in thermal contact such that
# 
# eqnarray169
# 
# and the total Hamiltonian is just  tex2html_wrap_inline620
# 
# Since system 2 is infinitely large compared to system 1, it acts as an infinite heat reservoir that keeps system 1 at a constant temperature T without gaining or losing an appreciable amount of heat, itself. Thus, system 1 is maintained at canonical conditions, N,V,T.
# 
# The full partition function  tex2html_wrap_inline626 for the combined system is the microcanonical partition function
# 
# displaymath173
# 
# Now, we define the distribution function,  tex2html_wrap_inline628 of the phase space variables of system 1 as
# 
# displaymath175
# 
# Taking the natural log of both sides, we have
# 
# displaymath177
# 
# Since  tex2html_wrap_inline630 , it follows that  tex2html_wrap_inline632 , and we may expand the above expression about  tex2html_wrap_inline634 . To linear order, the expression becomes
# 
# eqnarray179
# 
# where, in the last line, the differentiation with respect to  tex2html_wrap_inline636 is replaced by differentiation with respect to E. Note that
# 
# eqnarray184
# 
# where T is the common temperature of the two systems. Using these two facts, we obtain
# 
# eqnarray191
# 
# Thus, the distribution function of the canonical ensemble is
# 
# displaymath197
# 
# The prefactor  tex2html_wrap_inline642 is an irrelevant constant that can be disregarded as it will not affect any physical properties.
# 
# The normalization of the distribution function is the integral:
# 
# displaymath200
# 
# where Q(N,V,T) is the canonical partition function. It is convenient to define an inverse temperature  tex2html_wrap_inline646 . Q(N,V,T) is the canonical partition function. As in the microcanonical case, we add in the ad hoc quantum corrections to the classical result to give
# 
# displaymath204
# 
# The thermodynamic relations are thus,
# 
# Hemlholtz free energy:
# displaymath210
# 
# To see that this must be the definition of A(N,V,T), recall the definition of A:
# 
# displaymath213
# 
# But we saw that
# 
# displaymath215
# 
# Substituting this in gives
# 
# displaymath219
# 
# or, noting that
# 
# displaymath222
# 
# it follows that
# 
# displaymath229
# 
# This is a simple differential equation that can be solved for A. We will show that the solution is
# 
# displaymath232
# 
# Note that
# 
# displaymath235
# 
# Substituting in gives, therefore
# 
# displaymath241
# 
# so this form of A satisfies the differential equation.
# 
# Other thermodynamics follow:
# 
# Average energy:
# eqnarray244
# 
# Pressure:
# displaymath250
# 
# Entropy:
# eqnarray257
# 
# Heat capacity at constant volume:
# displaymath269
