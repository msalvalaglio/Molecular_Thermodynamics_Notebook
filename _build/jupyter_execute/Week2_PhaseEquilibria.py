#!/usr/bin/env python
# coding: utf-8

# # Liquid/Vapour Equilibrium conditions: deriving the Raoult Equation 
# 
# To express the partial pressure of component $i$ in a vapour phase in equilibrium with a multicomponent liquid phase we shall begin by expressing the fundamental condition that has to be fulfilled by all species at thermodynamic equilibrium at constant $T$ and $P$: 
# 
# $$
# \mu_i^L(T,P,\vec{x})=\mu_i^V(T,P,\vec{y})
# $$(PVAPeq:fund)
# 
# where $\mu_i^L(T,P,\vec{x})$ and $\mu_i^V(T,P,\vec{y})$ are the chemical potentials for specie $i$ in the liquid and vapour phase, respectively. 
# 
# For mixtures of ideal gases the chemical potential is conveniently expressed in differential form as a function of the partial pressure of component $i$, $p_i$ as $d\mu_i=RTd\ln{p_i}$. This result is obtained integrating the Gibbs-Duhem equation (_can you show it?_). 
# 
# In analogy with this expression in the case of real gases fugacity $f_i$ is introduced, defined as the partial pressure multiplied by a fugacity coefficient $\phi_i$, and leading to the differential expression $d\mu_i=RTd\ln{f_i}$. 
# 
# Integrating this expression between the phases in equilibrium yields: 
# 
# $$
# \int_{L,T,P,\vec{x}}^{V,T,P,\vec{y}}d\mu_i=\int_{L,T,P,\vec{x}}^{V,T,P,\vec{y}}RTd\ln{f_i}
# $$
# 
# $$
# \mu_i^V(T,P,\vec{y})-\mu_i^L(T,P,\vec{x})=RT\ln\frac{f_i^V(T,P,\vec{y}}{f_i^L(T,P,\vec{x})}
# $$(PVAPeq:chempot)
# 
# Hence, the equilibrium condition expressed by Eq. {eq}`PVAPeq:fund` can be conveniently stated as: 
# 
# $$
# f_i^V(T,P,\vec{y})=f_i^L(T,P,\vec{x})
# $$(PVAPeq:EQUI)
# 
# The left-hand term can be written as: 
# 
# $$
# f_i^V(T,P,\vec{y})=\phi(T,P,\vec{y})Py_i
# $$(PVAPeq:vapour)
# 
# where $\phi(T,P,\vec{y})$ is the fugacity coefficient, $P$ the pressure and $y_i$ the molar fraction in the vapour phase. 
# 
# The right hand side term can instead be written as: 
# 
# $$
# f_i^L(T,P,\vec{x})=\gamma_i(T,P,\vec{x})x_i\phi(T,P^o(T))P^o(T)e^{\frac{v_i(P-P^o(T))}{RT}}
# $$(PVAPeq:liquid)
# 
# where $x_i$ is the molar fraction of component $i$ in the liquid phase, $\gamma_i(T,P,\vec{x})$ is the activity coefficient for component i in the liquid mixture, $P^o(T)$ is the equilibrium vapour pressure of the pure component $i$, $\phi(T,P^o(T))$ is the fugacity coefficient for the pure component $i$ at $T$ and $P^o(T)$, and the term $e^{\frac{v_i(P-P^o(T))}{RT}}$ is the Poynting correction, which captures the difference in fugacity of the pure liquid component associated with the difference in pressure between $P$ and $P^o(T)$. 
# 
# 
# 
# Introducing Eq.{eq}`PVAPeq:vapour` and Eq. {eq}`PVAPeq:liquid` in Eq.{eq}`PVAPeq:EQUI`, while considering negligible the Poynting correction, and ideal gas approximation applicable yields: 
# 
# $$
# Py_i=\gamma_i(T,P,\vec{x})P^o(T)x_i
# $$(PVAPeq:final)
# 
# It should be noted that introducing the further hypothesis of ideal liquid mixture yields the Raoult law: 
# 
# $$
# Py_i=P^o(T)x_i
# $$(PVAPeq:Raoult)
# 
