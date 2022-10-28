#!/usr/bin/env python
# coding: utf-8

# # Week 4 Exercise: Second Virial Coefficient from Interatomic potentials. 
# 
# Using the expression: 
# 
# $$
# B_{ii}=2\pi{N_A}\int_0^\infty\left(1-e^{-\Gamma_{ii}/RT}\right)r^2dr
# $$
# 
# Compute the second virial coefficient for the following interatomic potentials: 
# 
# - $\Gamma_{ii}=0$ 
# 
# - $\Gamma_{ii}=0$ for $r> \sigma$ ; $\infty$ for $r \leq \sigma$
# 
# - $\Gamma_{ii}=0$ for $r>R\sigma$ ; $-\epsilon$ for $\sigma<r\leq{R}\sigma$ ; $\infty$ for $r \leq \sigma$
# 

# In[1]:


import matplotlib.pyplot as plt 
from matplotlib import cm
import numpy as np

figure=plt.figure()
axes = figure.add_axes([0.1,0.1,1.5,1.5])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
axes.set_xlabel('$\sigma$', fontsize=16);
axes.set_ylabel('B',fontsize=16);
axes.set_xlim([0.5,20]);
axes.set_ylim([-250,50]);

sigma=np.linspace(1E-4,1E-3,10)

NA=6E23

rho=1;

B=2/3*NA*np.pi*sigma**3

Z=1+B*rho

print(Z)
plt.scatter(sigma,Z)


# In[ ]:




