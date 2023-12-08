import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm 
import img
from scipy.ndimage import convolve
from scipy.optimize import curve_fit
import pandas as pd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Serif"],
})

N = int(input('Lattice size: '))

def random_lattice(N):
    '''
    Generate random initial condition.
    '''
    np.random.seed(10)
    return np.random.choice([1,0], [N,N])

def glider(i, j, lattice):
    '''
    Add a glider with top left cell at (i,j).
    '''

    glider = np.array([[0, 1, 0],
                      [0, 0, 1],
                      [1, 1, 1]])
    lattice[i:i+3, j:j+3] = glider

    return lattice

def blinker(i, j, lattice):
    '''
    Add a blinker with top left cell at (i,j).
    '''

    blinker = np.array([[0, 1, 0],
                       [0, 1, 0],
                       [0, 1, 0]])
    lattice[i:i+3, j:j+3] = blinker

    return lattice

def absorber(i, j, lattice):
    '''
    Add an absorber with top left cell at (i, j).
    '''

    absorber = np.array([[0, 0, 0],
                        [0, 1, 0],
                        [0, 0, 0]])
    lattice[i:i+3, j:j+3] = absorber

    return lattice

def beehive(i, j, lattice):
    '''
    Add a beehive with top left cell at (i, j).
    '''

    beehive = np.array([[0, 1, 0],
                       [1, 0, 1],
                       [1, 0, 1],
                       [0, 1, 0]])
    lattice[i-2:i+2, j-1:j+2] = beehive

    return lattice

def evolution(frames, img, lattice):
    '''
    Any live cell with less than 2 live neighbours dies.
    Any live cell with 2 or 3 live neighbours lives on to the next step.
    Any live cell with more than 3 live neighbours dies.
    Any dead cell with exactly 3 live neighbours comes to life.
    '''

    new_lattice = lattice.copy()

    # Convolve the current state of the lattice with a kernel
    # to sum up all living cells around the centre.
    padding = np.pad(lattice, pad_width=1, mode='wrap')
    kernel = np.array([[1,1,1], [1,0,1], [1,1,1]]) 
    live_lattice = convolve(padding, kernel, mode='constant')[1:-1,1:-1] 
    
    # Evolution conditions.
    new_lattice[np.where((new_lattice==1) & (live_lattice < 2))] = 0
    new_lattice[np.where((new_lattice==1) & (live_lattice > 3))] = 0
    new_lattice[np.where((new_lattice==0) & (live_lattice==3))] = 1

    if img != None:
        img.set_data(new_lattice)
    lattice[:] = new_lattice[:]

    return lattice
    
def simulation():
    print('Setting animation parameters...')
    # State can be random, glider, blinker, absorber, beehive.
    state = input('State: ')

    # Set animation parameters.
    step = int(input('Update interval: '))

    # Declare lattice.
    lattice = np.array([])

    # Set lattice for different initial conditions.
    if state == 'random':
        np.random.seed(None)
        lattice = random_lattice(N)

    elif state == 'glider':
        lattice = np.zeros((N,N))
        glider(1, 1, lattice)

    elif state == 'blinker':
        lattice = np.zeros((N,N))
        blinker(int(N/2), int(N/2), lattice)

    # Plot animation for the evolution of the system.
    fig, ax = plt.subplots()
    img = ax.imshow(lattice, interpolation='nearest')
    ani = animation.FuncAnimation(fig, evolution, fargs=(img, lattice, ),
                                  frames=10,
                                  interval=step,
                                  save_count=step)
    plt.show()

def equilibration():
    sims = 1000 # number of simulations
    steps = 10000 # number of steps
    time = [] # time list

    # Go through simulations, monitoring the active sites as a function of time.
    for i in tqdm(range(sims)):
        lattice = random_lattice(N)
        active_sites = []
        
        for step in range(steps):
            # Evolve the lattice and add active sites to list.
            new_lattice = evolution(frames=10, img=None, lattice=lattice)
            active_sites.append(np.sum(new_lattice))

            # Measure time when system has stopped evolving.
            if step > sims and np.all((active_sites[step-sims:step] == active_sites[step])) == True:
                time.append(step)
                break
    
    # Plot a histogram of the active sites over time.
    plt.hist(time, 50)
    plt.title('Distribution of Equilibration Times')
    plt.savefig('equilibration.png', dpi=1000)
    plt.show()

def centre_of_mass(masses):
    '''
    Calculate centre of mass from mass grid.
    '''
    return np.array((np.average(masses[0]), np.average(masses[1])))

def linear(m, x, b):
    '''
    Model linear function to fit to CM data.
    '''
    return m * x + b

def glider_cm():
    '''
    Track centre of mass for glider as a function of time.
    '''

    lattice = np.zeros((N,N))
    lattice = glider(1, 1, lattice)

    # Mass grid and initial centre of mass.
    masses = np.nonzero(lattice)
    cm = centre_of_mass(masses)

    # Lists for CM, time, and distance from (i,j).
    cm_list = []
    time_list = []
    distances = []

    time = 0
    # Track the CM as the lattice evolves over time.
    for i in range(N*2):
        time_list.append(time)
        time += 1
        
        new_lattice = evolution(frames=10, img=None, lattice=lattice)
        new_masses = np.nonzero(new_lattice)

        cm = centre_of_mass(new_masses)
        cm_list.append(cm.tolist())
        distances.append(np.linalg.norm(cm))

    data = [time_list, cm_list, distances]
    np.savetxt('./cm_data.txt', np.column_stack(data), header='time cm distance')
    
def plt_speed():
    glider_cm()

    df = np.loadtxt('cm_data', skiprows=1)
    x, y = np.arange(len(df[:,1])), df[:,1]

    # Calculate optimal values for CM data and covariance of these values.
    popt, pcov = curve_fit(linear, x, y)

    # Calculate standard deviation of CM data.
    err = np.sqrt(np.diag(pcov))

    # Plot the CM distance from (i,j) as a function of time.
    plt.errorbar(x, y, fmt='x')  
    plt.plot(x, linear(x, *popt), 'g')
    plt.title('CM distance from origin vs time')
    plt.xlabel('Time')
    plt.ylabel('CM distance')

    plt.savefig('cm_distance', dpi=1000)

if __name__ == '__main__':
    simulation()
    equilibration()
    plt_speed()