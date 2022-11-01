#!/usr/bin/env python
# coding: utf-8

# # From Intermolecular Interactions to Dynamics at the Atomistic Scale

# ## System Setup

# In[2]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt 
from matplotlib import cm
import numpy as np

## Initial conditions
# positions
x0=np.array([1, 2, 3, 3, 3, 2, 1, 1]);       
y0=np.array([1, 1, 1, 2, 3, 3, 3, 2]);

# timestep
dt=0.05;

#mass of the particles
m=np.array([1,1,1,1,1,1,1,1]);

#number of iterations
final_time=50;
NS=final_time/dt; 
nsteps=np.round(NS); 

#Interatomic potential constants
k=50.0; # Harmonic oscillator constant
req=1; # Harmonic oscillator equilibrium distance
HS=10; # Repulsive soft potential 

# Topology
M=np.array([[0, 1, 0, 0, 0, 0, 0, 0],
           [1, 0, 1, 0, 0, 0, 0, 0],
           [0, 1, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 1, 0, 0],
           [0, 0, 0, 0, 1, 0, 1, 0],
           [0, 0, 0, 0, 0, 1, 0, 1],
           [0, 0, 0, 0, 0, 0, 1, 0]]); 
k=k*M;

# Random initial velocities
v0=0.25*(np.random.rand(2,8)-0.5);


# ## Integration and Visualization

# ### Initialise the system

# In[3]:


get_ipython().run_cell_magic('capture', '', '%matplotlib inline\n\n# Setup figure for plotting the trajectory\n\nfigure, axes = plt.subplots(figsize=(5, 5))\nplt.xticks(fontsize=14)\nplt.yticks(fontsize=14)\naxes.set_xlim([-5,5]);\naxes.set_ylim([-5,5]);\n\n## Compute a trajectory with the Verlet Algorithm\n# Initialise positions at t-dt\nxp=x0;\nyp=y0;\n\n# Position at time t\nx=xp+v0[0,:]*dt;\ny=yp+v0[1,:]*dt;\n\n# Position at time t+dt\nxnew=np.zeros(np.shape(x0));\nynew=np.zeros(np.shape(x0));\n\n# time\ntime=np.arange(0,nsteps);\ncolor=iter(cm.gist_heat(np.linspace(0,1,np.size(time)+1)))\nxx=np.zeros((np.size(time),np.size(x)));xx[0]=x0\nyy=np.zeros((np.size(time),np.size(y)));yy[0]=y0\ntime[0]=0;\ntime[1]=time[0]+dt;\n\n# Initialise Energy Potential and Kinetic\nPOT=np.zeros(np.shape(time));\nKIN=np.zeros(np.shape(time));\n')


# ### Compute Trajectory

# In[4]:


# Compute trajectory
for timestep in np.arange(1,nsteps): #Cycle over timesteps
    c=next(color)                    #Update color for static representation
    timestep=int(timestep)           #Make sure timestep is an integer
    
    # Initialise force vectors
    fx=np.zeros(np.size(x0));  
    fy=np.zeros(np.size(x0)); 
    
    # Compute distances and the interparticle forces for every pair of particles
    for i in np.arange(0,np.size(x0)):
        for j in np.arange(i+1,np.size(x0)):
            
            r=np.sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])); #Distance    
               
            cx=-(k[i,j]*(r-req)-2*HS/(np.power(r,3)))*((x[i]-x[j]))/r;  #Pairwise force in x
            cy=-(k[i,j]*(r-req)-2*HS/(np.power(r,3)))*((y[i]-y[j]))/r;  #Pairwise force in y
                
            fx[i]=fx[i]+cx;      #update total x-component of the force on particle i
            fx[j]=fx[j]-cx;      #update total x-component of the force on particle j 

            fy[i]=fy[i]+cy;      #update total y-component of the force on particle i
            fy[j]=fy[j]-cy;      #update total y-component of the force on particle j 
           
       
    #Verlet integration
    for i in np.arange(0,np.size(x0)):
        xnew[i]=2*x[i]-xp[i]+(dt*dt)*fx[i]/m[i]; # new position (x-component)
        ynew[i]=2*y[i]-yp[i]+(dt*dt)*fy[i]/m[i]; # new position (y-component)
    
    # Compute velocity     
    vx=(xnew-xp)/2/dt;
    vy=(ynew-yp)/2/dt;
    v=np.sqrt(np.power(vx,2)+np.power(vy,2)); 
    
    # Reassign positions
    xp=x; yp=y; x=xnew-np.mean(xnew); y=ynew-np.mean(ynew);
    
    # Static representation of the trajectory
    line, = axes.plot(x,y,marker='o',color=c,markersize=10,linestyle='-')
    
    ## Store trajectory for animation 
    xx[timestep]=x;
    yy[timestep]=y;


# ### Visualization of the trajectory

# In[26]:


get_ipython().run_cell_magic('capture', '', "%matplotlib inline\nfrom matplotlib.animation import FuncAnimation\nfrom matplotlib import animation, rc\nfrom IPython.display import HTML\n\nfig, ax = plt.subplots(figsize=(8, 8))\nline, = ax.plot([]) \nax.set_xlim(-5, 5)\nax.set_ylim(-5, 5)\nline, = ax.plot([], [], lw=2, marker='o', markersize=45, markerfacecolor=(0.8, 1.0, 0.8, 0.5),\n             markeredgewidth=1,  markeredgecolor=(0, 0, 0, .5), linestyle='--',color='red')\n# initialization function: plot the background of each frame\ndef init():\n    line.set_data([], [])\n    return (line,)\n\ndef animate(frame_num):\n    x=xx[frame_num,:]\n    y=yy[frame_num,:]\n    line.set_data((x, y))\n    return (line,)\n\n# call the animator. blit=True means only re-draw the parts that have changed.\nanim = animation.FuncAnimation(fig, animate, init_func=init,\n                               frames=np.size(np.arange(1,nsteps)), interval=50);\n")


# In[27]:


HTML(anim.to_html5_video())


# In[ ]:




