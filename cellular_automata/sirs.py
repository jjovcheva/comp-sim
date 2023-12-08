import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import seaborn as sns
from tqdm import tqdm
import random
import pandas as pd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Serif"],
})

N = int(input('Lattice size: '))
n_sites = N ** 2 

# Assign values for susceptible, infected, and recovered/immune states.
S, I, R = 0, 1, -1

def random_lattice(N):
    '''
    Generate random initial condition.
    '''
    return np.random.choice([0,1,-1], [N,N])

def dynamics(lattice, p1, p2, p3):
    '''
    Update the system based on the dynamics of the SIRS model.
    '''
    # N = int(input('Lattice size: '))
    # p1 = float(input('p(S -> I): '))
    # p2 = float(input('p(I -> R): '))
    # p3 = float(input('p(R -> S): '))
    
    # Choose a random number between 0 and 1.
    r = np.random.rand()

    # Choose a lattice site (i, j) at random.
    i, j = random.choices(np.arange(N), k=2)

    # S -> I with p=p1 if at least 1 nearest neighbour of the site is I.
    if lattice[i, j] == S and p1 > r: 
        if (lattice[(i+1)%N, j] == I or\
            lattice[(i-1)%N, j] == I or\
            lattice[i, (j-1)%N] == I or\
            lattice[i, (j+1)%N] == I):
        
                lattice[i, j] = I

    # I -> R with probability p2.
    elif lattice[i,j] == I and p2 > r:
        lattice[i, j] = R

    # R -> S with probability p3.
    elif lattice[i, j] == R and p3 > r:
        lattice[i, j] = S

    return lattice

def simulation():
    # Choose evolution.
    evolution = input('System evolution: \n\
                      a: absorbing state \n\
                      b: dynamic equilibrium \n\
                      c: wave\n')
    
    if evolution == 'a': # absorbing state
        p1, p2, p3 = 0.5, 0.6, 0.1
    elif evolution == 'b': # dynamic equilibrium
        p1, p2, p3 = 0.5, 0.5, 0.5
    elif evolution == 'c': # wave
        p1, p2, p3 = 0.8, 0.1, 0.01
    
    steps = 1000
    # Generate random lattice.
    lattice = random_lattice(N)

    # Evolve lattice according to chosen dynamics.
    # Show animation.
    plt.figure
    im = plt.imshow(lattice, animated=True)
    for step in tqdm(range(steps)):
        for i in range(n_sites):
            lattice = dynamics(lattice, p1, p2, p3)
        
        plt.cla()
        plt.title(step)
        im = plt.imshow(lattice, cmap='magma', animated=True)
        plt.draw()
        plt.pause(0.0001)

def bootstrap(data, n):
    '''
    Measure bootstrap errors on data.
    '''
    res = np.random.choice(data, (n, 1000))
    f = np.var(res, axis=0)

    return np.var(f)**(1/2)

def phase():
    '''
    Calculate the average number of infected sites <I>. 
    Keeping p2 fixed at 0.5, track phase changes in the p1-p3 plane.

    :params p1_i, p3_i: initial values of p1, p3
    :params p1_f, p3_f: final values of p1, p3
    :param n: number of samples in p1, p3
    '''
    p2 = 0.5
    steps = 1101 # allowing 100 steps as equilibration time

    # Create array for p1 and p3 with n evenly spaced intervals each.
    p1_array = np.linspace(0., 1., 21)
    p3_array = np.linspace(0., 1., 21)

    # Combine p1, p3 to create probability grid.
    p_grid = [(p1, p3) for p1 in p1_array for p3 in p3_array]

    # Lists for average number of infected sites, variance, and error.
    I_avg = []
    I_var = []

    # Set up lattice and probability grid.
    for i in tqdm(range(len(p_grid))):
        I_sites = [] # list of infected sites
        p1, p3 = p_grid[i]
        lattice = random_lattice(N)

        # Evolve system and track number of infected sites over time.
        for step in tqdm(range(steps)):
            for i in range(n_sites):
                lattice = dynamics(lattice, p1, p2, p3)
            if step > 100:
                I_sites.append(np.count_nonzero(lattice==I))
  
        # Update lists.
        I_avg.append(np.mean(I_sites))
        I_var.append(np.var(I_sites))

    I_avg, I_var = np.array(I_avg)/n_sites, np.array(I_var)

    # Save data to text.
    data = [I_avg, I_var]
    np.savetxt('I_data.txt', np.column_stack(data), header='I_avg I_var')
    np.savetxt('p1-p3.txt', p_grid, header='p1 p3')

    # Build 2d grid from 1d p1 and p3 arrays.
    XX, YY = np.meshgrid(p1_array, p3_array)

    # Plot phase diagram of p1-p3 plane.
    fig = plt.figure()
    ax = plt.contourf(XX, YY, I_avg.reshape((len(p1_array), len(p3_array))), cmap='magma')
    plt.colorbar()
    plt.savefig('phase_plot', dpi=1000)

    # Plot contours of variance of the number of infected sites.
    fig = plt.figure()
    ax = plt.contourf(XX, YY, I_var.reshape((len(p1_array),len(p3_array))), cmap='magma')
    plt.colorbar()
    plt.savefig('var_plot', dpi=1000)

    return I_avg, I_var, p1_array, p3_array

def phase_cut():
    '''
    Calculate the average number of infected sites <I>. 
    Keeping p2 fixed at 0.5, track phase changes in the p1-p3 plane.

    :params p1_i, p3_i: initial values of p1, p3
    :params p1_f, p3_f: final values of p1, p3
    :param n: number of samples in p1, p3

    '''
    p2 = 0.5
    steps = 10101 # allowing 100 steps as equilibration time

    # Create array for p1 and p3 with n evenly spaced intervals each.
    p1_array = np.linspace(0.2, 0.5, 21)
    p3_array = np.linspace(0.5, 0.5, 1)

    # Combine p1, p3 to create probability grid.
    p_grid = [(p1, p3) for p1 in p1_array for p3 in p3_array]

    # Lists for average number of infected sites, variance, and error.
    I_avg = []
    I_var = []
    I_errs = []

    # Set up lattice and probability grid.
    for i in tqdm(range(len(p_grid))):
        I_sites = [] # list of infected sites
        p1, p3 = p_grid[i]
        lattice = random_lattice(N)

        # Evolve system and track number of infected sites over time.
        for step in tqdm(range(steps)):
            for i in range(n_sites):
                lattice = dynamics(lattice, p1, p2, p3)
            if step > 100:
                I_sites.append(np.count_nonzero(lattice==I))

        # Estimate bootstrap errors on data where p1 or p3 fixed.
        if (len(p1_array) == 1 or len(p3_array) == 1):
            I_errs.append(bootstrap(np.array(I_sites), 100))
  
        # Update lists.
        I_avg.append(np.mean(I_sites))
        I_var.append(np.var(I_sites))

    I_avg, I_var, I_errs = np.array(I_avg)/n_sites, np.array(I_var), np.array(I_errs)/n_sites

    # Save data to text.
    data = [I_avg, I_var, I_errs]
    np.savetxt('I_data_cut.txt', np.column_stack(data), header='I_avg I_var I_errs')

    # Plot variance with error bars.
    fig = plt.figure()
    ax = plt.plot(p1_array, I_var)
    ax = plt.errorbar(p1_array, I_var, yerr=I_errs, fmt='x', ecolor='darkmagenta', color='m')
    plt.xlabel('$p_1$')
    plt.ylabel('$I_\mathrm{var}$')
    plt.savefig('cut_var_plot', dpi=1000)

def f_Im():
    '''
    Model permanent immunity, tracking the average infected fraction.
    '''

    # Create array for fraction of immune agents.
    f_array = np.arange(0., 1.1, 0.1)

    # Lists for average infected fraction and errors.
    I_avg = []

    # Loop over array, turn a certain number of sites immune, and count infected sites.
    for im_frac in tqdm(f_array):
        I_single = []

        for i in tqdm(range(5)):
            I_sites = []

            lattice = random_lattice(N)
            i = np.random.randint(0, lattice.shape[0], int(lattice.size*im_frac))
            j = np.random.randint(0, lattice.shape[0], int(lattice.size*im_frac))
            lattice[i, j] = 3.

            steps = 1101

            # Evolve lattice and track number of infected sites over time.
            for step in tqdm(range(steps)):
                for site in range(n_sites):
                    lattice = dynamics(lattice, 0.5, 0.5, 0.5)

                if step > 100:
                    I_sites.append(np.count_nonzero(lattice==I))
                
            I_single.append(np.mean(I_sites))
 
        I_avg.append(np.mean(I_single))

    I_avg = np.array(I_avg)/n_sites
 
    np.savetxt('f_Im_data.txt', I_avg, header='I_avg')

    fig = plt.figure()
    ax = plt.plot(f_array, I_avg)
    ax = plt.errorbar(f_array, I_avg, fmt='x', ecolor='darkmagenta', color='m')
    plt.xlabel('$f_\mathrm{Im}$')
    plt.ylabel('$I_\mathrm{avg}$')
    plt.savefig('f_Im_plot', dpi=1000)

if __name__ == '__main__':
    #simulation()
    #phase()
    phase_cut()
    #f_Im()