The Cahn-Hilliard calculation should be run by typing ``python3 cahn-hilliard.py``
The only input is the initial value of $\phi$ (included here are plots for $\phi_0=0$ and $\phi_0=0.5$)
An animation will come up which is updated every 1000 steps, for a total of 100,000 steps.

The Poisson calculation should be run by typing ``python3 poisson.py``.
The inputs are the field ('electric' or 'magnetic') and the mode of calculation ('jacobi', 'gaussian' or 'sor').
The optimal value of omega for the electric case was found to be 1.9.
