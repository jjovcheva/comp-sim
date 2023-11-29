import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm 
import sys

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Serif"],
})

# Global variables
dx = 1.
e0 = 1.
N = 50
steps = 2000
tol = 1e-3
field = input('Field: ')

def pot_calc(potential, rho):
    '''
    Function to calculate potential.
    '''
    return (1/6) * (np.roll(potential, 1, axis=0) +\
                      np.roll(potential, -1, axis=0) +\
                      np.roll(potential, 1, axis=1) +\
                      np.roll(potential, -1, axis=1) +\
                      np.roll(potential, 1, axis=2) +\
                      np.roll(potential, -1, axis=2) + rho)

def checkerboard(shape):
    '''
    Create checkerboard, dividing lattice into black and white nodes.
    '''
    mask =  np.bool_(np.indices(shape).sum(axis=0) % 2)
    return mask, ~mask

def parameters(field):
    '''
    Choose parameters for boundary value problem.
    :param: field: 'electric' or 'magnetic'

    :return: shape: shape of lattice
    :return: s: slice of lattice
    :return: varphi: potential
    :return: reset: reset potential
    :return: rho: charge density
    :return: padding: padding for boundary conditions
    :return: active_site: active region for boundary conditions
    :return: black, white: checkerboard pattern
    '''
    shape = np.full(3, N+2)
    s = slice(1, N+1)

    varphi, reset, rho = np.zeros((3, *shape))
    black, white = checkerboard(shape)

    if field == 'electric':
        rho[N//2, N//2, N//2] = 1
        padding = ((1,1), (1,1), (1,1))
        active_site = (s, s, s)
    
    elif field == 'magnetic':
        rho[N//2, N//2] = 1
        padding = ((1,1), (1,1), (0,0))
        active_site = (s, s, )

    return shape, s, varphi, reset, rho, padding, active_site, black, white

shape, s, varphi, reset, rho, padding, active_site, black, white = parameters(field)

def calc_field(field):
    '''
    Calculate field from potential.
    :return: f_x, f_y, f_z: field components
    '''

    f_x = (np.roll(varphi, 1, axis=0) - np.roll(varphi, -1, axis=0))[s, s, N//2].ravel()/(-2 * dx)
    f_y = (np.roll(varphi, 1, axis=1) - np.roll(varphi, -1, axis=1))[s, s, N//2].ravel()/(-2 * dx)
    f_z = (np.roll(varphi, 1, axis=2) - np.roll(varphi, -1, axis=2))[s, s, N//2].ravel()/(-2 * dx)

    return f_x, f_y, f_z

def jacobi():
    global varphi
    '''
    Solve boundary value problem using Jacobi algorithm.
    :return: distance between successive estimates of varphi
    '''

    varphi_0 = np.copy(varphi)
    varphi = pot_calc(varphi, rho)
    
    varphi = np.pad(varphi[s,s,s], pad_width=1)
    return np.sum(np.abs(varphi - varphi_0))

def gauss_seidel():
    '''
    Solve boundary value problem using Gauss-Seidel algorithm.
    :return: err: distance between successive estimates of varphi
    '''
    global varphi

    varphi_0 = np.copy(varphi)
    varphi = pot_calc(varphi, rho)
    
    # Reverse the black updates.
    varphi[black] = varphi_0[black]

    # Boundary conditions.
    varphi = np.pad(varphi[active_site], padding, mode='constant', constant_values=0)
    white_updates = np.copy(varphi)

    # Calculate potential.
    varphi = pot_calc(varphi, rho)

    # Reverse the white updates.
    varphi[white] = white_updates[white]
    varphi_f = np.pad(varphi[active_site], padding, mode='constant', constant_values=0)

    err = np.sum(np.abs(varphi_f - varphi_0))

    return err

def sor(omega=1):
    '''
    Solve boundary value problem using SOR algorithm.
    :param: omega: relaxation parameter
    :return: err: distance between successive estimates of varphi
    '''

    global varphi

    varphi_0 = np.copy(varphi)
    varphi = pot_calc(varphi, rho)
    
    # Overrelax.
    varphi *= omega
    varphi += (1 - omega) * varphi_0

    # Reverse the black updates.
    varphi[black] = varphi_0[black]

    # Boundary conditions.
    varphi = np.pad(varphi[active_site], padding, mode='constant', constant_values=0)
    white_updates = np.copy(varphi)

    varphi = pot_calc(varphi, rho)
    
    # Overrelax.
    varphi *= omega
    varphi += (1 - omega) * varphi_0

    # Reverse the white updates.
    varphi[white] = white_updates[white]

    # Boundary conditions.
    varphi_f = np.pad(varphi[active_site], padding, mode='constant', constant_values=0)
    err = np.sum(np.abs(varphi_f - varphi_0))

    return err

def save_data(field):

    pos = np.array(np.meshgrid(np.arange(1, N+1), np.arange(1, N+1))).T.reshape(-1,2)
    dis = np.linalg.norm(pos - N//2, axis=1)
    potentials = varphi[1:N+1, 1:N+1, N//2].ravel()

    if field == 'electric':
        E_x, E_y, E_z = -1 * np.array(calc_field(field))
        np.savetxt('electric_data.txt', np.c_[pos, dis, potentials, E_x, E_y, E_z], header='pos dis potentials E_x E_y E_z')

        # Plot electric field.
        E_mag = np.sqrt(E_x**2 + E_y**2 + E_z**2)
        plt.yscale('linear')
        plt.quiver(pos.T[0], pos.T[1], E_x/E_mag, E_y/E_mag, E_mag, angles='xy', scale_units='xy', scale=1)
        plt.title('Electric Field')
        plt.savefig('e_field.png', dpi=800)

        plt.clf()
        plt.scatter(dis, E_mag, marker='x', s=15, c='g')
        plt.xlabel('Distance')
        plt.ylabel('$E_\mathrm{mag}$')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Electric field strength vs. distance')
        plt.savefig('e_mag.png', dpi=800)

    elif field == 'magnetic':
        B_y, B_x, _= calc_field(field)
        np.savetxt('magnetic_data.txt', np.c_[pos, dis, potentials, B_x, -B_y], header='pos dis potentials B_x -B_y')

        # Plot magnetic field.
        B_mag = np.sqrt(B_x**2 + B_y**2)
        plt.yscale('linear')
        plt.quiver(pos.T[0], pos.T[1], B_x/B_mag, -B_y/B_mag, B_mag, angles='xy', scale_units='xy', scale=1)
        plt.title('Magnetic Field')
        plt.savefig('b_field.png', dpi=800)

        plt.clf()
        plt.scatter(dis, B_mag, marker='x', s=15, c='b')
        plt.xlabel('Distance')
        plt.ylabel('$B_\mathrm{mag}$')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Magnetic field strength vs. distance')
        plt.savefig('b_mag.png', dpi=800)

    np.savetxt('field_2d.txt', varphi[0:N+1, 0:N+1, N//2])

def plot_slice(mode):
    # Plot slice of the potential.
    plt.imshow(varphi[1:N+1, N//2, 1:N+1], interpolation='gaussian', cmap='magma')
    plt.colorbar()
    plt.savefig('slice_%s.png' %mode, dpi=800)

mode = input('Mode: ')
if mode == 'jacobi':
    # Calculate using Jacobi method and plot slice of the field.

    for sweep in tqdm(range(steps)):
        if np.isclose(jacobi(), 0, atol=tol):
            break
    plot_slice(mode)

elif mode == 'gaussian':
    # Calculate using Gauss-Seidel method, plot slice of the field,
    # and save data.

    for sweep in tqdm(range(steps)):
        err = gauss_seidel()
        if np.isclose(err, 0, atol=tol):
            break
    plot_slice(mode)
    save_data(field)

elif mode == 'sor':
    # Calculate using SOR method, plot sweeps vs. omega, and find
    # optimal value of omega.

    iters = 50
    omegas = np.linspace(1, 2, iters, endpoint=False)
    data = []

    for i in tqdm(range(iters)):
        varphi = np.copy(reset)
        # print(omegas[i])
        for sweep in tqdm(range(steps)):
            err = sor(omega=omegas[i])
            if np.isclose(err, 0, atol=0.1):
                data += [[omegas[i], sweep]]
                break

    data = np.array(data)
    
    # Find and print the optimal value of omega.
    condition = np.where(data[:,1] == np.min(data[:,1]), data[:,0], 0)
    print('The optimal value of omega is: ', condition[condition != 0])

    # Plot sweeps vs. omega.
    plt.clf()
    plt.plot(data[:,0], data[:,1], color='darkmagenta') 
    plt.xlabel('$\omega$')
    plt.ylabel('Sweeps')
    plt.title('Sweeps vs. $\omega$ for %s field' %field)
    plt.savefig('sor_%s.png' %field, dpi=800)
    np.savetxt('sor_data_%s.txt' %field, data, header='omega sweeps')
        

