# 2D-Ice-Shell-Creep
Viscous creep in a spherical shell for power-law rheology

Step 1: Adjust constants in "Constants"
=

Step 2: Select creep mechanism and basal viscosity (eta_0) in "Rheology"
=

creep= [1], [2], or [3]:

[1]= newtonian diffusion creep (n=1)

[2]= grain boundary sliding (n=1.8)

[3]= dislocation creep (n=4)

10^12 < eta_0 <10^14 recommended

Step 3: Set dimensionless grid spacing and timestep in "Mesh"
=

suggested timesteps for eta0=10^12,Ntheta=51,Nz=101: tconv/100 (n=4), 5tconv (n=1), 2tconv (n=1.8) 

set trans_rel to estimated percentage of viscous deformation. 0.25-0.5 is recommended to start

*tconv is 1 yr in seconds

Step 4: Prescribe initial condition in "initial thickness perturbation":
=

set Hmin to minimum thickness, 
Hpole to initial polar/maximum thickness,
and Heq to predicted equilibrium thickness

set dtheta to half the width of thinned region, in radians 

*dtheta of pi/2 corresponds to global-scale variation

Step 5: Run "Simulation"
=
