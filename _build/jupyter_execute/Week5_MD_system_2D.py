#!/usr/bin/env python
# coding: utf-8

# # MD in a system of pseudo-molecules in 2D
# 
# ## System Setup

# In[1]:


## Useful functions
verlet=lambda r, r_past, force, mass, dt:  2*r-r_past+(dt**2)*force/mass
forcebox=lambda x, boxx,boxk: np.greater(np.abs(x),boxx)*(-boxk)*x


#Define the system's box
boxx=10 #x dimension of the simulation' box
boxy=10 #y dimension of the simulation' box
boxk=1  #k constant for harmonic repulsive force

#Number of particles
N=8

#mass of the particles
m=np.ones(N)

# Topology
M=np.array([[0, 1, 0, 0, 0, 0, 0, 1],
           [1, 0, 1, 0, 0, 0, 0, 0],
           [0, 1, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 1, 0, 0],
           [0, 0, 0, 0, 1, 0, 1, 0],
           [0, 0, 0, 0, 0, 1, 0, 1],
           [1, 0, 0, 0, 0, 0, 1, 0]]); 

M_gas_diatomic=np.array([[0, 1, 0, 0, 0, 0, 0, 0],
           [1, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 1, 0, 0],
           [0, 0, 0, 0, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 1],
           [0, 0, 0, 0, 0, 0, 1, 0]]); 




######## "Force Field" Parameters #######
HS=1; # Repulsive soft potential 
k=25.0; # Harmonic oscillator constant
req=1; # Harmonic oscillator equilibrium distance
KAPPA=k*M_gas_diatomic
epsilon=0;


##Use this function to implement different potentials
def forceij(xi,xj,yi,yj,HS,KAPPA,req,epsilon): 
        r=np.sqrt((xi-xj)**2+(yi-yj)**2); #Distance      
                
        #Ideal Gas
        dVdr=0 
        #Repulsive Wall 
        # dVdr=-12*HS/(np.power(r,13))
        #Repulsive Wall + Harmonic potential
        #... ... ...
        #Lennard Jones Potential
        #... .... .... ... ....         
        cx=-dVdr*((xi-xj))/r;  #Pairwise component of the force in x
        cy=-dVdr*((yi-yj))/r;  #Pairwise component of the force in y
             
        return [cx,cy]


def print_progress(iteration, total, bar_length=50):
    progress = (iteration / total)
    arrow = '*' * int(round(bar_length * progress))
    spaces = ' ' * (bar_length - len(arrow))
    print(f'\r|{arrow}{spaces}| {int(progress * 100)}% | ', end='', flush=True)


# In[ ]:


## Set the initial Conditions
# Random initial positions
x0=(np.random.rand(N)*2*boxx)-(boxx);     #Initial position in x
y0=(np.random.rand(N)*2*boxy)-(boxy);     #Initial position in y 

# Random initial velocities
v0=(np.random.rand(2,N)-0.5); # Initial random velocitites

## Define the timestep and the total time
dt=0.005; # Timestep
total_time=100;  # Total simulation time
nsteps=int(total_time/dt); # Total number of steps

## Initialise vectors 
time=np.zeros(nsteps)


# ## Integration and Visualization

# ### Initialise the system

# In[ ]:


## Compute a trajectory with the Verlet Algorithm
# Initialise positions at t-dt
xp=x0;
yp=y0;

# Position at time t
x=xp+v0[0,:]*dt;
y=yp+v0[1,:]*dt;

# Position at time t+dt
xnew=np.zeros(N);
ynew=np.zeros(N);

# time
time=np.arange(0,nsteps);
time[0]=0;
time[1]=time[0]+dt;

## Initialize verctors for plotting 
xx=np.zeros((np.size(time),N));xx[0]=x0
yy=np.zeros((np.size(time),N));yy[0]=y0


# ### Compute Trajectory

# In[ ]:


## |------------------|
## |Compute trajectory|
## |------------------|
for timestep in np.arange(1,nsteps): #Cycle over timesteps
    timestep=int(timestep)           #Make sure timestep is an integer
    
    # Initialise force vectors
    fx=np.zeros(N);  
    fy=np.zeros(N); 
    
    # Cycle over all particles
    for i in np.arange(0,N):
        fx[i]+=forcebox(x[i],boxx,boxk)
        fy[i]+=forcebox(y[i],boxy,boxk)
        for j in np.arange(i+1,N):
            
            [cx,cy]=forceij(x[i],x[j],y[i],y[j],HS,KAPPA,req,epsilon)
            
            fx[i]=fx[i]+cx;      #update total x-component of the force on particle i
            fx[j]=fx[j]-cx;      #update total x-component of the force on particle j 

            fy[i]=fy[i]+cy;      #update total y-component of the force on particle i
            fy[j]=fy[j]-cy;      #update total y-component of the force on particle j 

        xnew[i]=verlet(x[i],xp[i],fx[i],m[i],dt) # new position (x-component)
        ynew[i]=verlet(y[i],yp[i],fy[i],m[i],dt); # new position (y-component)

    print_progress(timestep,nsteps)  

    # Reassign positions
    xp=x; yp=y; x=xnew+1-1; y=ynew+1-1;

    ## Store trajectory for animation 
    xx[timestep]=x;
    yy[timestep]=y;


# ### Visualization of the trajectory

# In[ ]:


get_ipython().run_cell_magic('capture', '', "## Display the trajectory\n%matplotlib inline\nfrom matplotlib.animation import FuncAnimation\nfrom matplotlib import animation, rc\nfrom IPython.display import HTML\n\nfig, ax = plt.subplots(figsize=(10, 10))\nline, = ax.plot([]) \nax.set_xlim(-boxx, boxx)\nax.set_ylim(-boxy, boxy)\nline, = ax.plot([], [], lw=2, marker='o', markersize=40, markerfacecolor=(0.8, 1.0, 0.8, 0.5),\n             markeredgewidth=1,  markeredgecolor=(0, 0, 0, .5), linestyle=' ',color='red')\n# initialization function: plot the background of each frame\ndef init():\n    line.set_data([], [])\n    return (line,)\n\ndef animate(frame_num):\n    x=xx[frame_num,:]\n    y=yy[frame_num,:]\n    line.set_data((x, y))\n    return (line,)\n\n# call the animator. blit=True means only re-draw the parts that have changed.\nanim = animation.FuncAnimation(fig, animate, init_func=init,\n                               frames=np.arange(1,int(nsteps),50), interval=50);\n")


# In[ ]:


HTML(anim.to_jshtml())


# In[ ]:


HTML(anim.to_jshtml())

