"""
CMod Project B: velocity Verlet time integration of
a molecular system of particles interacting via
the Lennard-Jones potential.

Takes in an input file with defined number of particles,
density, temperature, timestep, numstep, and cutoff radius.

Outputs 4 files: 
An XYZ file with particle trajectories
A file with the system's total, kinetic, and potential
energies over time
A file with the mean squared displacement over time
A file with the average radial distribution

The potential is U = 4 * ((1/r12 ^ 12)-(1/r12 ^ 6)).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import mdutilities
import pbc
from p3d_a import Particle3D as p3d

def lj_force_pot(separation, cutoff):
    
    """
    Method to return the force and potential
    energy between 2 particles interacting via
    the Lennard-Jones potential.

    The force is given by
    F = 48((1/r12^14)-(1/(2 * r12^8)))*(r1-r2)

    The potential energy is given by
    U = 4 * ((1/r12 ^ 12)-(1/r12 ^ 6))

    :param separation: separation between 2
    particles in the system
    :param cutoff: cutoff radius past which
    the force between particles is 0
    :return: force acting on particle, potential 
    energy between a particle pair
    """

    r12_v = separation
    r12 = np.linalg.norm(separation)
    if r12 < cutoff:
        force = 48 * ((1 / r12 ** 14)-(1 / (2 * r12 ** 8))) * r12_v
        potential = 4 * ((1/r12 ** (12))-(1/r12 ** 6))
    else: 
        force = 0
        potential = 4 * ((1/cutoff** (12))-(1/cutoff ** 6))
    return force, potential

def net_force_pot(separation, cutoff):
    """
    Method to return the force and potential
    energies between all particles in a list. 
    Takes in the separation and for N particles
    adds the force acting on them from the 
    lj_force_pot function. Applies the same
    action for the potential energy.

    :param separation: separation between 2
    particles in the system
    :param cutoff: cutoff radius past which
    the force between particles is 0
    :return force list, potential: force array, 
    potential energy of the system
    """

    N = len(separation)
    force_list = np.zeros((N,3))
    potential = 0
    for i in range(N):
        for j in range(N):
            if i != j:
                for_pot = lj_force_pot(separation[i][j], cutoff)
                force_list[i] += for_pot[0]
                potential += for_pot[1]

    return force_list, potential/2

# Main method
def main():
    # Read name of input file from command line
    if len(sys.argv)!=6:
         print("Wrong number of arguments.")
         quit()
    else:
        infile_name = sys.argv[1]
        outfile1_name = sys.argv[2]
        outfile2_name = sys.argv[3]
        outfile3_name = sys.argv[4]
        outfile4_name = sys.argv[5]

    # Open 4 output files
    # Trajectory output in XYZ format
    outfile1 = open(outfile1_name +".xyz", 'w+')
    # Output file with total, kinetic + potential energies
    outfile2 = open(outfile2_name, 'w+') 
    # Output files with MSD and RDF
    outfile3 = open(outfile3_name, 'w+')
    outfile4 = open(outfile4_name, 'w+')

    # Open input file
    file = open(infile_name + '.txt','r')
    f = file.readlines()
    parameter = f[1].split()

    # Take in user-defined parameters from input file
    n = int(parameter[0]) # Number of particles in the system
    rho = float(parameter[1]) # Density
    T = float(parameter[2]) # Temperature
    dt = float(parameter[3])
    numstep = int(parameter[4])
    cutoff = float(parameter[5]) # Cutoff radius
    time = 0.0
    dr = 0.1 
    
    # Initialise Particle3D instances for n particles
    p3d_list = []
    for i in range(n):
        p3d_list.append(p3d.new_particle(f[0],i))

    # Close input file
    file.close()

    # Use mdutilities to set box size and initial positions/velocities
    box_size, none = mdutilities.set_initial_positions(rho, p3d_list)
    mdutilities.set_initial_velocities(T, p3d_list)

    # Get initial particle separation, force, and potential energy
    separation, sep_list = p3d.p3d_sep(p3d_list, box_size)
    force, sys_pe = net_force_pot(separation, cutoff)

    # Calculate MSD
    ref_list = [p.pos for p in p3d_list]
    msd = p3d.msd(p3d_list, ref_list, box_size)
    
    # Bin initial separations into a histogram
    his, bin_edges = p3d.rdf(separation, box_size, dr)

    # Get initial particle energy
    sys_ke = p3d.sys_kinetic(p3d_list)
    energy = p3d.sys_kinetic(p3d_list) + sys_pe

    # Initialise data lists for plotting later
    time_list = [time]
    energy_list = [energy]
    ke_list = [sys_ke]
    pe_list = [sys_pe]
    msd_list = [msd]
    sep_list = [sep_list]

    # Start the time integration loop
    for i in range(numstep):
        # Update particle positions
        p3d.update_pos_2nd_list(p3d_list, dt, force, box_size)
        
        # Update force
        sep_new, sep_mod = p3d.p3d_sep(p3d_list, box_size)
        force_new, sys_pe = net_force_pot(sep_new, cutoff)

        # Compute MSD
        msd = p3d.msd(p3d_list, ref_list, box_size)

        # Bin separations into a histogram
        his, bin_edges = p3d.rdf(separation, box_size, dr)

        # Update particle velocities by averaging current and new forces
        p3d.update_vel_list(p3d_list, dt, 0.5*(force+force_new))
        
        # Re-define force value
        force = force_new

        # Increase time
        time += dt

        # Write out particle trajectories in XYZ format
        
        outfile1.write("{0}\n".format(n))
        outfile1.write("Point = {}\n".format(i+1))
        for p in p3d_list:
            outfile1.write(str(p)+"\n") 
                
        # Write out particle energy over time
        energy = p3d.sys_kinetic(p3d_list) + sys_pe
        outfile2.write("{0:5.3f} {1:12.8f}  {2:12.8f}  {3:12.8f}\n".format(time, energy, p3d.sys_kinetic(p3d_list), sys_pe)) 

        # Write out mean squared displacement over time
        outfile3.write("{0:1.3f} {1:12.8f} \n".format(time, msd))     

        # Append information to data lists
        time_list.append(time)
        energy_list.append(energy)
        ke_list.append(p3d.sys_kinetic(p3d_list))
        pe_list.append(sys_pe)
        msd_list.append(msd)
        sep_list.append(sep_mod)

        # Histogram update
        his, bin_edges = p3d.rdf(sep_list, box_size, dr)
        # Expected RDF for homogeneous system
        rho_0 = 4 * np.pi * rho * (np.square(bin_edges[:-1])) * dr
        # Calculate radial distribution function
        rdf = his / (rho_0 * n * numstep)

    # Write out radial distribution function by separation
    for i,j in zip(bin_edges[:-1], rdf):
        outfile4.write("{0:f} {1:f} \n".format(i,j))

    # Post-simulation:
    # Close output files
    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()     


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

