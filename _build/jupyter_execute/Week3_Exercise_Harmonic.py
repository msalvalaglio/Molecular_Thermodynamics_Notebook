#!/usr/bin/env python
# coding: utf-8

# # Week 3 Exercise: The harmonic oscillator. 
# 
# Given a 1D Harmonic oscillator with characteristic frequency $\omega=5$Hz, mass of 20, starting at the initial position of 1, with an initial momentum of 1. 
# 
# - Can You define the coordinates of the Phase space of the Harmonic oscillator? 
# 
# - How would you represent the motion of the Harmonic oscillator in Phase Space? 
# 
# - What is the most probable state for the Harmonic oscillator? 

# In[8]:


get_ipython().run_line_magic('matplotlib', 'notebook')
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt 
import numpy as np

## Parameters
m=1; #mass
k=2; # Harmonic Constant
v0=1; 
x0=1; 

om=np.sqrt(k/m)

t=np.linspace(0,10,100)

print(t)

x=x0*np.cos(np.prod(om,t))+v0*sin(om*t)

print(x)

plt.plot(t,x); plt.show()


# In[ ]:




