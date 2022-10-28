#!/usr/bin/env python
# coding: utf-8

# # The Compressibility Factor
# 
# The phase behavior of fluids at moderate to high pressure requires the definition of an equation of state that captures deviations from the ideal gas behaviour. 
# 
# To this aim we introduce a compressibility factor: 
# 
# $$
# Z=\frac{Pv}{RT}
# $$
# 
# for an ideal gas $Z=1$, $P=RT/v$, and $v=RT/P$. 
# 
# # Fugacity and compressibility factor
# 
# As discussed when working on Vapour/Liquid Equilibria, the chemical potential of a real fluid can be computed from its fugacity $f_i$. 
# The ratio between the fugacity of a pure component and its pressure can be written as: 
# 
# $$
# RT\ln\left(\frac{f_i}{P}\right)_{pure}=RT\ln{\phi_{i,pure}}=\int_0^P\left(v_i-\frac{RT}{P}dP\right) 
# $$
# 
# this can be recast as a function of the compressibility factor as: 
# 
# $$
# \ln{\phi_{i,pure}}=\int_0^P\left(\frac{Z_i-1}{P}dP\right) 
# $$
# 
# hence being able to caputre the behavior of the compressibility factor in a real fluid opens up the possibility of characterising equilibria beyond the ideal gas approximation. 
# 

# ## The Virial equation of state
# 
# The virial equation of state is based on the expression of the compressibility factor as a power series of the density $\rho=v^{-1}$. 
# 
# $$
# Z=1+\rho{B}+\rho^2{C}+\rho^3{D}+...
# $$
# 
# This expansion is usually truncated at the third term. Each term, represents the contribution to the compressibility factor that comes from interactions between two, three, four, etc bodies. 
# For pure components the virial coefficients are independent from pressure or density and only function of the system's temperature
# 
# This equation also provides a direct point of contact between molecular and macroscopic properties. In fact, experimentally the virial expansion coefficients can be obtained by estimating the following limits: 
# 
# $$
# B=\lim_{\rho\rightarrow{0}}\left(\frac{\partial{Z}}{\partial{\rho}}\right)
# $$
# 
# $$
# C=\lim_{\rho\rightarrow{0}}\frac{1}{2}\left(\frac{\partial^2{Z}}{\partial{\rho^2}}\right)
# $$
# 
# practically this can be done by estimating as the interceipt and slope of the following equation: 
# 
# $$
# v\left(\frac{Pv}{RT}-1\right)=B+\frac{C}{v}
# $$
# 
# 
