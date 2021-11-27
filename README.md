# NavierStokes_Solver_MPI_OpenMP
Efficient hybrid schemes to solve Navier Stokes equation with periodic boundary conditions using pseudo spectral methods.

The programs numerically solves the 3D incompressible Navier-Stokes on a Cubic Domain [0,2pi]x[0,2pi]x[0,2pi] using pseudo-spectral methods and Implicit Midpoint rule timestepping. The numerical solution is compared to an exact solution reported by Shapiro 

Analytical Solution:

u(x,y,z,t)=-0.25*(cos(x)sin(y)sin(z)+sin(x)cos(y)cos(z))exp(-t/Re)

v(x,y,z,t)= 0.25*(sin(x)cos(y)sin(z)-cos(x)sin(y)cos(z))exp(-t/Re)

w(x,y,z,t)= 0.5*cos(x)cos(y)sin(z)exp(-t/Re)
