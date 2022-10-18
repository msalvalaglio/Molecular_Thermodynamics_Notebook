#!/usr/bin/env python
# coding: utf-8

# I. THE MICROSCOPIC LAWS OF MOTION
# Consider a system of N classical particles. The particles a confined to a particular region of space by a “container” of volume V . The particles have a finite kinetic energy and are therefore in constant motion, driven by the forces they exert on each other (and any external forces which may be present). At a given instant in time t, the Cartesian positions of the particles are r1 (t), ..., rN (t). The time evolution of the positions of the particles is then given by Newton’s second law of motion:
# mi ̈ri =Fi(r1,...,rN)
# where F1,...,FN are the forces on each of the N particles due to all the other particles in the system. The notation
#  ̈ri = d2ri/dt2.
# N Newton’s equations of motion constitute a set of 3N coupled second order differential equations. In order to solve these, it is necessary to specify a set of appropriate initial conditions on the coordinates and their first time derivaties, {r1(0), ..., rN (0), r ̇1(0), ..., r ̇N (0)}. Then, the solution of Newton’s equations gives the complete set of coordinates and velocities for all time t.

# In[ ]:


The idea of ensemble averaging can also be expressed in terms of an average over all such microstates (which comprise the ensemble). A given macroscopic property, A, and its microscopic function a = a(x), which is a function of the positions and momenta of a system, i.e. the phase space vector, are related by

1 􏰂N
A = ⟨a⟩ensemble = N a(xλ)
λ=1 where xλ is the microstate of the λth member of the ensemble.

However, recall the original problem of determining the microscopic detailed motion of each individual particle in a system. In reality, measurements are made only on a single system and all the microscopic detailed motion is present. However, what one observes is still an average, but it is an average over time of the detailed motion, an average that also washes out the microscopic details. Thus, the time average and the ensemble average should be equivalent, i.e.

1􏰃T
A = ⟨a⟩ensemble = lim dt a(x(t))
This statement is known as the ergodic hypothesis. A system that is ergodic is one for which, given an infinite amount of time, it will visit all possible microscopic states available to it (for Hamiltonian dynamics, this means it will visit all points on the constant energy hypersurface). No one has yet been able to prove that a particular system is truly ergodic, hence the above statement cannot be more than a supposition. However, it states that if a system is ergodic, then the ensemble average of a property A(x) can be equated to a time average of the property over an ergodic tra jectory.

This statement is known as the ergodic hypothesis. A system that is ergodic is one for which, given an infinite amount of time, it will visit all possible microscopic states available to it (for Hamiltonian dynamics, this means it will visit all points on the constant energy hypersurface). No one has yet been able to prove that a particular system is truly ergodic, hence the above statement cannot be more than a supposition. However, it states that if a system is ergodic, then the ensemble average of a property A(x) can be equated to a time average of the property over an ergodic tra jectory.


