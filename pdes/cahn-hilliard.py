import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm 
from scipy.ndimage import laplace # can use this instead of rolling function, it makes little difference.

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Serif"],
})

# Global variables
N = 100
M = 0.1
a = 0.1
kappa = 0.1
dx = 1.
dt = 1.

phi_0 = float(input('Enter initial value of phi: '))

def rolling(var):
    '''
    Rolling function to calculate the Laplacian of a variable.
    '''
    return (np.roll(var, 1, axis=0) +
         np.roll(var, -1, axis=0) +
         np.roll(var, 1, axis=1) +
         np.roll(var, -1, axis=1) -
         np.multiply(var, 4))

def ch(phi):
    '''
    Update lattice according to Cahn-Hilliard equation.

    :param: phi: compositional order parameter
    :return: dphi: change in compositional order parameter
    '''
    mu = -a * phi + a * np.power(phi, 3) - (kappa/dx ** 2) * rolling(phi) 
    
    dphi = M * (dt/dx ** 2) * rolling(mu)
    
    return dphi

def f(phi):
    '''
    Calculate free energy density of the system.

    :param: phi: compositional order parameter
    :return: f: free energy density
    '''
    dphi_0 = np.gradient(phi, axis=0) **2
    dphi_1 = np.gradient(phi, axis=1) **2

    f = np.sum(-(a/2)*np.power(phi, 2) + (a/4)*np.power(phi, 4) + kappa/2*(dphi_0 + dphi_1))
    
    return f

def animate():
    '''
    Animate the Cahn-Hilliard equation on a lattice.    
    '''

    # Initialise lattice with some random noise.
    phi = np.random.uniform(phi_0 + 0.1, phi_0 - 0.1, (N, N))

    # Number of timesteps.
    steps = int(1e+5)+1

    # Initialise figure.
    plt.figure()
    plt.imshow(phi, vmin=-1, vmax=1, cmap='magma')
    plt.colorbar()

    energies = np.empty((steps//100)+1)

    for i in tqdm(range(steps)):
        # Update phi.
        phi += ch(phi)

        if i%1000 == 0:
            # Update figure every 1000 steps.
            plt.cla()   
            plt.title(i)
            plt.imshow(phi, vmin=-1, vmax=1, animated=True, interpolation='gaussian', cmap='magma')
            plt.draw()
            plt.pause(0.0001)
        
        if i%100 == 0:
            # Save free energy density every 100 steps.
            energies[i//100] = f(phi)

    step = np.arange(0, int(1e+5)+100, 100)
    np.savetxt('f_%s.txt' %phi_0, np.c_[step, energies], header='timestep f')

def plot_energy():
    '''
    Plot free energy density vs timestep.
    '''

    plt.clf()
    x, y = np.loadtxt('f_%s.txt' %phi_0, unpack=True)
    plt.plot(x, y, color='darkmagenta')
    plt.xlabel('Timestep')
    plt.ylabel('Free energy density')
    plt.title('Free energy density vs timestep for $\phi_0$ = %s' %phi_0)
    plt.savefig('f_%s.png' %phi_0, dpi=800)

animate()
#plot_energy()

