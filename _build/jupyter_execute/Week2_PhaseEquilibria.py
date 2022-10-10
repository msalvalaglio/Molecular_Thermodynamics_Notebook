#!/usr/bin/env python
# coding: utf-8

# # Topic 1: Liquid/Vapour Equilibrium conditions: deriving the Raoult Equation 
# 
# To express the partial pressure of component $i$ in a vapour phase in equilibrium with a multicomponent liquid phase we shall begin by expressing the fundamental condition that has to be fulfilled by all species at thermodynamic equilibrium at constant $T$ and $P$: 
# 
# $$
# \mu_i^L(T,P,\vec{x})=\mu_i^V(T,P,\vec{y})
# $$(PVAPeq:fund)
# 
# where $\mu_i^L(T,P,\vec{x})$ and $\mu_i^V(T,P,\vec{y})$ are the chemical potentials for specie $i$ in the liquid and vapour phase, respectively. 
# 
# For mixtures of ideal gases the chemical potential is conveniently expressed in differential form as a function of the partial pressure of component $i$, $p_i$ as $d\mu_i=RTd\ln{p_i}$.  
# In analogy with this expression in the case of real gases fugacity $f_i$ is introduced, leading to the differential expression $d\mu_i=RTd\ln{f_i}$. 
# Integrating this expression between the phases in equilibrium yields: 
# 
# $$
# \int_{L,T,P,\vec{x}}^{V,T,P,\vec{y}}d\mu_i=\int_{L,T,P,\vec{x}}^{V,T,P,\vec{y}}RTd\ln{f_i}
# $$
# $$
# \mu_i^V(T,P,\vec{y})-\mu_i^L(T,P,\vec{x})=RT\ln\frac{f_i^V(T,P,\vec{y}}{f_i^L(T,P,\vec{x})}
# $$(PVAPeq:chempot)
# 
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
# where $x_i$ is the molar fraction of component $i$ in the liquid phase, $\gamma_i(T,P,\vec{x})$ is the activity coefficient for component i in the liquid mixture, $P^o(T)$ is the equilibrium vapour pressure of the pure component $i$, $\phi(T,P^o(T))$ is the fugacity coefficient for the pure component $i$ at $T$ and $P^o(T)$, and the term $e^{\frac{v_i(P-P^o(T))}{RT}}$ is the Poynting correction. For a discussion on the origin of the Poynting correction please refer to the notes on osmosis. 
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

# # Topic 2: Liquid/Vapour Equilibrium conditions: deriving the Raoult Equation 
# 
# Let us begin by writing the solid-liquid equilibrium conditions at temperature $T$ and pressure $P$, between a pure crystalline phase and a multicomponent liquid phase characterised by composition $\vec{x}$. 
# Equilibrium conditions is stated by the equality of chemical potentials: 
# 
# $$
# \mu_i^s(T,P)=\mu_i^\ell(T,P,\vec{x})
# $$(SOLeq1)
# 
# which is equivalent to the equality of the fugacity in the two phases: 
# 
# $$
# f_i^s(T,P)=f_i^\ell(T,P,\vec{x})
# $$(SOLeq2)
# 
# Focusing on the right hand of Eq.{eq}`SOLeq2`, the fugacity in the liquid phase can be rewritten as a function of the fugacity of the pure liquid phase as follows:
# 
# $$
# f_i^\ell(T,P,\vec{x})=f_i^\ell(T,P)x_i\gamma_i(T,P,\vec{x})
# $$(SOLeq3)
# 
# where $x_i$ is the molar fraction of the solute and $\gamma_i(T,P,\vec{x})$ is the activity coefficient of that specie in the liquid phase.  
# 
# Let us now focus on the left hand side of Eq. {eq}`SOLeq2`. The fugacity of the pure solid phase can also be written as a function of the fugacity of the pure liquid: 
# 
# $$
# \int^{T,P,s}_{T,P,\ell}d\mu_i=\int^{T,P,s}_{T,P,\ell}RTd\ln{f_i}
# $$(SOLeq4)
# 
# $$
# \Delta{\mu}_{i,\ell\rightarrow{s}}=RTln\left(\frac{f_i^s(T,P)}{f_i^\ell(T,P)}\right)
# $$(SOLeq5)
# 
# $$
# f_i^s(T,P)=f_i^\ell(T,P)e^{\frac{\Delta{\mu}_{i,\ell\rightarrow{s}}}{RT}}
# $$(SOLeq6)
# 
# The change in chemical potential can be written for a pure substance as the partial molar change in free energy $\Delta{g_i(T,P)}$:
# 
# $$
# \Delta\mu_{i,\ell\rightarrow{s}}=\Delta{g}_i^s(T,P)=\Delta{h}_i(T,P)-T\Delta{s}_i(T,P)
# $$(SOLeq7)
# 
# 
# The molar enthalpy change is: 
# 
# $$
# \Delta{h}_i(T,P)=-\Delta{h_{fus}\left(T_f\right)}+\int_T^{T_f}Cp^{\ell}dT+\int_T^{T_f}Cp^{s}dT\simeq-\Delta{h_{fus}\left(T_f\right)}+\left(Cp^{\ell}-Cp^{s}\right)\left(T_f-T\right)
# $$(SOLeq8)
# 
# while the molar entropy change is written as: 
# 
# $$
# \Delta{s}_i(T,P)=-\frac{\Delta{h_{fus}\left(T_f\right)}}{T_f}+\int_T^{T_f}\frac{Cp^{\ell}}{T}dT+\int_T^{T_f}\frac{Cp^{s}}{T}dT \simeq -\frac{\Delta{h_{fus}}}{T_f}+\left(Cp^{\ell}-Cp^{s}\right)\ln{\left(\frac{T_f}{T}\right)}
# $$(SOLeq9)
# 
# Defining $\Delta{Cp}_{\ell\rightarrow{s}}=\left(Cp^{\ell}-Cp^{s}\right)$
# 
# $$
# \Delta{g}_i^s(T,P)=-\Delta{h_{fus}\left(1-\frac{T}{T_f}\right)}+\Delta{Cp}_{\ell\rightarrow{s}}\left(T_f-T-T\ln{\left(\frac{T_f}{T}\right)}\right)
# $$(SOLeq10)
# 
# Since $\Delta{Cp}_{\ell\rightarrow{s}}$ is usually small compared the above expression is typically simplified to: 
# 
# $$
# \Delta{g}_i^s(T,P)\simeq-\Delta{h_{fus}\left(1-\frac{T}{T_f}\right)}
# $$(SOLeq11)
# 
# Inserting this equation into the equality of fugacities yields:
# 
# $$
# f_i^\ell(T,P)x_i\gamma_i(T,P,\vec{x})=f_i^\ell(T,P)e^{\frac{\Delta{\mu}_{i,\ell\rightarrow{s}}}{RT}}
# $$(SOLeq12)
# 
# hence the solubility of component $i$, $x_i$ is:  
# 
# $$
# x_i=\frac{1}{\gamma_i(T,P,\vec{x})}e^{-\frac{\Delta{h_{fus}}}{R}\left(\frac{1}{T}-\frac{1}{T_f}\right)}
# $$(SOLeq13)

# In[ ]:




