## Simulation of Argon atoms interacting via the Lennard-Jones potential ##
This project simulates Argon atoms under the Lennard-Jones potential at different phases. It was written in 2021 (my second year of university), so it is certainly not optimised, but it does what it needs to. 

The atoms' trajectory is propagated using Verlet time integration inside a simulation box that follows periodic boundary conditions and the minimum image convention.

The script `lj.py` takes one input file and produces 4 output files.
All values are in reduced units.

The code should be run by typing: `lj.py <input_handle>`
where `<input_handle>` is the name of the input txt file, e.g. 'gas'

Input - plain text file in the following format: 

    Line 1: `<label> <mass> <x> <y> <z> <vx> <vy> <vz>`
    
    Line 2: `<no. of particles> <density> <temperature> <timestep> <numstep> <cutoff radius>`

Output - four files with the following data:

    `traj.xyz`: particle trajectories in XYZ format
    `energy`: time and particle total, kinetic, and potential energies
    `msd`: time and mean squared displacement
    `rdf`: r and average radial distribution function
