#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
@author: s1850975
'''

import numpy as np
import pbc

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    p3d_sep - computes particle pair separations from a P3D list
    update_pos_2nd_list - second order position update for P3D list
    update_vel_list - velocity update for P3D list
    sys_kinetic - computes total K.E. of a P3D list
    com_velocity - computes total mass and CoM velocity of a P3D list
    msd - computes mean squared displacement of the system
    rdf - computes radial distribution function of the system
    """
    
    def __init__(self, label, mass, pos, vel):
        """
        Initialises a particle in 3D space.

        :param label: string w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """        
        self.label = str(label)
        self.mass = float(mass)
        self.pos = np.array(pos,float)
        self.vel = np.array(vel,float) 
      
    @staticmethod
    def new_particle(input_string, n):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one per particle in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param input_file: Readable file handle in the above format

        return Particle3D instance
        """
    
        tokens = input_string.split()
        label = tokens[0] + str(n)
        mass = tokens[1]
        pos = np.array([tokens[2],tokens[3], tokens[4]])
        vel = np.array([tokens[5], tokens[6], tokens[7]])
           
        return Particle3D(label, mass, pos, vel)

    @staticmethod
    def p3d_sep(p3d_list, box_size):
        """
        Calculates the separations for a list of P3D instances,
        following PBC and respecting MIC.
        Calculates the separation moduli.

        :param p3d_list: list in which each item is a P3D instance
        :param box_size: box size set by mdutilities

        :return separations, separation_mod: separations, separation moduli
        """
        
        N = len(p3d_list)
        separations = np.zeros((N,N,3))
        separation_mod = np.zeros((N,N))
        for i in range(N):
            for j in range(i):
                separations[i,j] = pbc.mic((p3d_list[i].pos - p3d_list[j].pos), box_size)
                separations[j,i] = -pbc.mic((p3d_list[i].pos - p3d_list[j].pos), box_size)
                separation_mod[i,j] = np.linalg.norm(pbc.mic((p3d_list[i].pos - p3d_list[j].pos), box_size))
                separation_mod[j,i] = -np.linalg.norm(pbc.mic((p3d_list[i].pos - p3d_list[j].pos), box_size))
        return separations, separation_mod

    def __str__(self): 
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """              
        xyz_string = str(self.label)+" "+str(self.pos[0])+" "+str(self.pos[1])+" "+str(self.pos[2])
        
        return xyz_string
        
    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m · vel^2
        """
        ke = 0.5 * self.mass * (np.linalg.norm(self.vel) ** 2)
        
        return ke 
    
    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m · vel
        """
        p = self.mass * self.vel
        
        return p
    
    def update_pos_1st(self, dt):
        """
        1st order position update

        :param dt: given timestep
        """
        self.pos = self.pos + dt * self.vel
    
    def update_pos_2nd(self, dt, f):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.pos = self.pos + (dt * self.vel) + ((dt ** 2) * f / (2 * self.mass))

    @staticmethod
    def update_pos_2nd_list(p3d_list, dt, f_list, box_size):
        """
        2nd order position update

        :param dt: timestep
        :param p3d_list: list in which each item is a P3D instance
        :param f_list: list of force instances
        """ 
        for p3d,f in zip(p3d_list,f_list):
            p3d.pos = pbc.pbc((p3d.pos + p3d.vel*dt + 0.5*(1/p3d.mass)*f*(dt**2)), box_size)
        
    def update_vel(self, dt, f):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        self.vel = self.vel + (dt * f / self.mass)

    @staticmethod
    def update_vel_list(p3d_list, dt, f_list):
        """
        Velocity list update

        :param p3d_list: list in which each item is a P3D instance
        :param dt: timestep
        :param f_list: force list
        """
        for p3d,f in zip(p3d_list,f_list):
            p3d.vel = p3d.vel + (dt * f/p3d.mass)
        
    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Returns the kinetic energy of the whole system

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: sum(1/2 m_i · vel_i^2)
        """
        sys_ke = 0
        
        for object in p3d_list:
            sys_ke += object.kinetic_e()
            
        return sys_ke
    
    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance

        :return total_mass: the total mass of the system 
        :return com_vel: centre-of-mass velocity
        """   
        total_mass = []
        total_p = []
        
        for object in p3d_list:
            total_mass.append(object.mass)   
            total_p.append(object.momentum())
      
        com_vel = np.sum(total_p) / np.sum(total_mass)
        
        return np.sum(total_mass), com_vel

    @staticmethod
    def msd(p3d_list, o_list, box_size):
        """
        Computes the mean squared displacement of the system.

        :param p3d_list: list in which each item is a P3D instance
        :param o_list: reference list with initial particle positions
        :param box_size: box size set by mdutilities

        :return d/N: mean squared displacement
        """
        N = len(p3d_list)
        d = 0
        for i, j in zip(p3d_list, o_list):
            d += np.square(np.linalg.norm(pbc.mic((i.pos-j), box_size)))
        
        return d / N

    @staticmethod
    def rdf(sep_list, box_size, dr):
        """
        Computes the radial distribution function of the system by
        binning pair separations into a histogram. For bins, uses
        values between 0.1 and the radius of the sphere that
        circumscribes the box defined by mdutilities.

        :param sep_list: list of particle pair separations
        :param box_size: box size set by mdutilities
        :param dr: bin width

        :return his, bin_edges: histogram and bin_edges
        """
        bins = np.arange(dr, (3**0.5/2)*box_size[0], dr)
        his, bin_edges = np.histogram(sep_list, bins = bins, density=False)

        return his, bin_edges

        
