This project explores partial differential equations on a lattice 
-- specifically, the Cahn-Hilliard (an example of an initial value problem)
and Poisson (example of a boundary value problem) equations.

## Cahn-Hilliard equation ##
The Cahn-Hilliard equation describes phase separation in a physical system,
such as water and oil emulsion. The equation is characterised by a
compositional order parameter, $\phi(\vec x, t)$ which depends on the position
$\vec x$ and time $t$, which will be positive if there is locally more
water than oil, and negative otherwise. The overall integral of $\phi$ must 
be conserved as the amount of each substance remains constant.

The order parameter obeys the Cahn-Hilliard equation
$$ \frac{\partial \phi(\vec{x}, t)}{\partial t} = M\nabla^2\mu(\vec{x}, t), $$

where $M$ is a positive constant and the chemical potential $\mu(\vec x, t)$
is $$\mu(\vec x, t) = -a\phi(\vec x, t)+a\phi(\vec x, t)^3 -\kappa\nabla^2\phi(\vec x, t),$$
where $a,\kappa$ are positive constants.

The free energy density $f$ associated with the order parameter is given by
$$f=-\frac{a}{2}\phi^2 +\frac{a}{4}\phi^4 +\frac{\kappa}{2}(\nabla\phi)^2.$$
Over time, we expect the system to evolve in a way that minimises the free energy density.

The Cahn-Hilliard calculation should be run by typing `python3 cahn-hilliard.py`.
The only input is the initial value of $\phi$ (included here are plots for $\phi_0=0$ and $\phi_0=0.5$). The animation is updated every 1000 steps, for a total of 100,000 steps.

## Poisson equation ##
The electrostatic potential $\varphi(\vec r)$ due to a scalar field of electric
charges $\rho(\vec r)$ is the solution to Poisson's equations
$$\nabla^2\varphi=-\frac\rho{\epsilon_0}$$
where $\epsilon_0$ is the vacuum dielectric constant. The electric field is given by $\vec E=-\nabla\varphi$. 

The evolution of the system is explored through the Jacobi, Gauss-Seidel, and
successive over-relaxation (SOR) algorithms.  
The Poisson calculation should be run by typing `python3 poisson.py`.
The inputs are the field ('electric' or 'magnetic') and the mode of calculation ('jacobi', 'gaussian' or 'sor').
The optimal value of $\omega$ for the electric case was found to be 1.9.
